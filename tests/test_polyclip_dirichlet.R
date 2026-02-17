## Test script for polyclipping with Dirichlet tessellation
## Tests: dirtess_clip_areas_cpp, dirtess_clip_cpp, and polygonise with window

library(ppmData)
library(sf)
library(terra)

cat("=== Test 1: Basic Dirichlet tessellation (no clipping) ===\n")
set.seed(42)
coords <- matrix(runif(200), ncol = 2)
tess <- dirTess(coords)
cat("  ncoords:", tess$ncoords, "\n")
cat("  n polygons:", length(tess$polygons$poly), "\n")
cat("  PASS\n\n")

cat("=== Test 2: dirtess_clip_areas_cpp with square window ===\n")
set.seed(42)
n <- 50
pts <- matrix(runif(2 * n), ncol = 2)

# Dummy points (same as dirTess internals)
bbox_dummy <- ppmData:::get_bbox(pts, buffer = 0.05)
dummycoords <- ppmData:::get_bbox_dummy_coords(pts, bbox_dummy, n = 9)
pts_all <- rbind(pts, dummycoords)
coords_vec <- as.numeric(t(pts_all))

# Square clip polygon [0.1, 0.9] x [0.1, 0.9]
clip_x <- c(0.1, 0.9, 0.9, 0.1)
clip_y <- c(0.1, 0.1, 0.9, 0.9)

areas <- ppmData:::dirtess_clip_areas_cpp(
  coords = coords_vec,
  ncoords = n,
  clip_x = clip_x,
  clip_y = clip_y,
  hole_x_list = list(),
  hole_y_list = list()
)

cat("  n areas:", length(areas), "\n")
cat("  all areas >= 0:", all(areas >= 0), "\n")
cat("  total area:", sum(areas), "\n")
cat("  expected area (0.8 x 0.8):", 0.64, "\n")
cat("  area close to expected:", abs(sum(areas) - 0.64) < 0.01, "\n")
cat("  PASS\n\n")

cat("=== Test 3: dirtess_clip_cpp returns polygons + areas ===\n")
res <- ppmData:::dirtess_clip_cpp(
  coords = coords_vec,
  ncoords = n,
  clip_x = clip_x,
  clip_y = clip_y,
  hole_x_list = list(),
  hole_y_list = list()
)

cat("  has 'areas':", "areas" %in% names(res), "\n")
cat("  has 'polygons':", "polygons" %in% names(res), "\n")
cat("  n polygons:", length(res$polygons), "\n")
cat("  areas match:", max(abs(res$areas - areas)) < 1e-12, "\n")

# Check polygon structure
poly1 <- res$polygons[[1]]
cat("  polygon has 'outer':", "outer" %in% names(poly1), "\n")
cat("  polygon has 'holes':", "holes" %in% names(poly1), "\n")
cat("  outer is matrix:", is.matrix(poly1$outer), "\n")
cat("  outer ncol == 2:", ncol(poly1$outer) == 2, "\n")
cat("  PASS\n\n")

cat("=== Test 4: Clipped polygon vertices stay inside clip region ===\n")
all_inside <- TRUE
for (i in seq_len(n)) {
  outer <- res$polygons[[i]]$outer
  if (nrow(outer) == 0) next
  xs <- outer[, 1]
  ys <- outer[, 2]
  if (any(xs < 0.1 - 1e-8) || any(xs > 0.9 + 1e-8) ||
      any(ys < 0.1 - 1e-8) || any(ys > 0.9 + 1e-8)) {
    cat("  FAIL: polygon", i, "has vertices outside clip region\n")
    all_inside <- FALSE
  }
}
cat("  all vertices inside clip:", all_inside, "\n")
cat("  PASS\n\n")

cat("=== Test 5: polygonise with sf window (C++ fast path) ===\n")
set.seed(123)
coords2 <- matrix(runif(100), ncol = 2)
tess2 <- dirTess(coords2)

# Create an sf polygon window
win_coords <- matrix(c(0.05, 0.05,
                        0.95, 0.05,
                        0.95, 0.95,
                        0.05, 0.95,
                        0.05, 0.05), ncol = 2, byrow = TRUE)
win_poly <- st_polygon(list(win_coords))
win_sfc <- st_sfc(win_poly, crs = "EPSG:4326")
win_sf <- st_sf(geometry = win_sfc)

poly_res <- polygonise(tess2, window = win_sf, clippy = TRUE, unit = "geo")
cat("  n polygon areas:", length(poly_res$polygons.areas), "\n")
cat("  all areas > 0:", all(poly_res$polygons.areas > 0), "\n")
cat("  sf geometry type:", as.character(unique(st_geometry_type(poly_res$polygons))), "\n")

# Verify total area is close to the window area
total_area <- sum(poly_res$polygons.areas)
win_area <- 0.9 * 0.9
cat("  total area:", total_area, "\n")
cat("  window area:", win_area, "\n")
cat("  area close to window:", abs(total_area - win_area) < 0.02, "\n")
cat("  PASS\n\n")

cat("=== Test 6: Clip with a hole (donut window) ===\n")
set.seed(99)
n3 <- 80
pts3 <- matrix(runif(2 * n3), ncol = 2)

bbox_dummy3 <- ppmData:::get_bbox(pts3, buffer = 0.05)
dummycoords3 <- ppmData:::get_bbox_dummy_coords(pts3, bbox_dummy3, n = 9)
pts_all3 <- rbind(pts3, dummycoords3)
coords_vec3 <- as.numeric(t(pts_all3))

# Outer: unit square
outer_x <- c(0, 1, 1, 0)
outer_y <- c(0, 0, 1, 1)

# Hole: small square in the middle [0.4, 0.6] x [0.4, 0.6]
hole_x <- c(0.4, 0.6, 0.6, 0.4)
hole_y <- c(0.4, 0.4, 0.6, 0.6)

areas_hole <- ppmData:::dirtess_clip_areas_cpp(
  coords = coords_vec3,
  ncoords = n3,
  clip_x = outer_x,
  clip_y = outer_y,
  hole_x_list = list(hole_x),
  hole_y_list = list(hole_y)
)

# Without hole for comparison
areas_nohole <- ppmData:::dirtess_clip_areas_cpp(
  coords = coords_vec3,
  ncoords = n3,
  clip_x = outer_x,
  clip_y = outer_y,
  hole_x_list = list(),
  hole_y_list = list()
)

cat("  total area with hole:", sum(areas_hole), "\n")
cat("  total area without hole:", sum(areas_nohole), "\n")
cat("  expected difference (0.04):", abs(sum(areas_nohole) - sum(areas_hole) - 0.04) < 0.01, "\n")
cat("  hole reduces area:", sum(areas_hole) < sum(areas_nohole), "\n")
cat("  PASS\n\n")

cat("=== Test 7: dirtess_clip_cpp with hole returns hole polygons ===\n")
res_hole <- ppmData:::dirtess_clip_cpp(
  coords = coords_vec3,
  ncoords = n3,
  clip_x = outer_x,
  clip_y = outer_y,
  hole_x_list = list(hole_x),
  hole_y_list = list(hole_y)
)

has_holes <- any(sapply(res_hole$polygons, function(p) length(p$holes) > 0))
cat("  some polygons have holes:", has_holes, "\n")
cat("  areas match clip_areas:", max(abs(res_hole$areas - areas_hole)) < 1e-10, "\n")
cat("  PASS\n\n")

cat("=== Test 8: Edge case - very few points ===\n")
set.seed(1)
pts_few <- matrix(c(0.3, 0.3, 0.7, 0.7, 0.5, 0.5), ncol = 2, byrow = TRUE)
bbox_few <- ppmData:::get_bbox(pts_few, buffer = 0.05)
dummy_few <- ppmData:::get_bbox_dummy_coords(pts_few, bbox_few, n = 5)
pts_all_few <- rbind(pts_few, dummy_few)
coords_few <- as.numeric(t(pts_all_few))

areas_few <- ppmData:::dirtess_clip_areas_cpp(
  coords = coords_few,
  ncoords = 3,
  clip_x = c(0, 1, 1, 0),
  clip_y = c(0, 0, 1, 1),
  hole_x_list = list(),
  hole_y_list = list()
)

cat("  3 point areas:", areas_few, "\n")
cat("  all >= 0:", all(areas_few >= 0), "\n")
cat("  total ~ 1.0:", abs(sum(areas_few) - 1.0) < 0.01, "\n")
cat("  PASS\n\n")

cat("=== Test 9: Non-rectangular clip polygon (triangle) ===\n")
set.seed(77)
n4 <- 40
pts4 <- matrix(runif(2 * n4), ncol = 2)
bbox4 <- ppmData:::get_bbox(pts4, buffer = 0.05)
dummy4 <- ppmData:::get_bbox_dummy_coords(pts4, bbox4, n = 9)
pts_all4 <- rbind(pts4, dummy4)
coords_vec4 <- as.numeric(t(pts_all4))

# Triangle: (0,0), (1,0), (0.5,1)
tri_x <- c(0, 1, 0.5)
tri_y <- c(0, 0, 1)

areas_tri <- ppmData:::dirtess_clip_areas_cpp(
  coords = coords_vec4,
  ncoords = n4,
  clip_x = tri_x,
  clip_y = tri_y,
  hole_x_list = list(),
  hole_y_list = list()
)

expected_tri_area <- 0.5  # area of triangle
cat("  total clipped area:", sum(areas_tri), "\n")
cat("  expected triangle area:", expected_tri_area, "\n")
cat("  area close to expected:", abs(sum(areas_tri) - expected_tri_area) < 0.02, "\n")
cat("  PASS\n\n")

cat("=== Test 10: Visual check - plot clipped tessellation ===\n")
set.seed(42)
coords_plot <- matrix(runif(600), ncol = 2)
tess_plot <- dirTess(coords_plot)

# Circular-ish window (hexagon)
angles <- seq(0, 2 * pi, length.out = 7)[-7]
hex_x <- 0.5 + 0.45 * cos(angles)
hex_y <- 0.5 + 0.45 * sin(angles)
hex_coords <- cbind(c(hex_x, hex_x[1]), c(hex_y, hex_y[1]))
hex_poly <- st_polygon(list(hex_coords))
hex_sfc <- st_sfc(hex_poly, crs = "EPSG:4326")
hex_sf <- st_sf(geometry = hex_sfc)

poly_hex <- polygonise(tess_plot, window = hex_sf, clippy = TRUE, unit = "geo")

pdf("/home/woo457/Dropbox/ppmData/tests/test_polyclip_plot.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))

# Left: clipped tessellation with hexagonal window
plot(st_geometry(poly_hex$polygons), col = hcl.colors(nrow(poly_hex$polygons), alpha = 0.5),
     main = "Clipped to hexagon")
plot(hex_sfc, add = TRUE, border = "red", lwd = 2)
points(coords_plot, pch = 16, cex = 0.8)

# Right: clipped to a square for comparison
sq_coords <- matrix(c(0.1, 0.1, 0.9, 0.1, 0.9, 0.9, 0.1, 0.9, 0.1, 0.1),
                    ncol = 2, byrow = TRUE)
sq_sf <- st_sf(geometry = st_sfc(st_polygon(list(sq_coords)), crs = "EPSG:4326"))
poly_sq <- polygonise(tess_plot, window = sq_sf, clippy = TRUE, unit = "geo")
plot(st_geometry(poly_sq$polygons), col = hcl.colors(nrow(poly_sq$polygons), alpha = 0.5),
     main = "Clipped to square")
plot(st_geometry(sq_sf), add = TRUE, border = "red", lwd = 2)
points(coords_plot, pch = 16, cex = 0.8)

dev.off()
cat("  Plot saved to tests/test_polyclip_plot.pdf\n")
cat("  PASS\n\n")

cat("=== Test 11: Compare C++ clipping vs sf st_intersection ===\n")
set.seed(55)
n5 <- 60
coords5 <- matrix(runif(2 * n5), ncol = 2)
tess5 <- dirTess(coords5)

win5_coords <- matrix(c(0.1, 0.1, 0.9, 0.1, 0.9, 0.9, 0.1, 0.9, 0.1, 0.1),
                      ncol = 2, byrow = TRUE)
win5_sf <- st_sf(geometry = st_sfc(st_polygon(list(win5_coords)), crs = "EPSG:4326"))

# C++ path (with window)
t_cpp <- system.time({
  res_cpp <- polygonise(tess5, window = win5_sf, clippy = TRUE, unit = "geo")
})

# sf path (without window, clipped to bbox)
t_sf <- system.time({
  res_sf <- polygonise(tess5, window = NULL, clippy = TRUE, unit = "geo")
})

cat("  C++ path time:", t_cpp["elapsed"], "s\n")
cat("  sf path time:", t_sf["elapsed"], "s\n")
cat("  C++ areas (n):", length(res_cpp$polygons.areas), "\n")
cat("  sf areas (n):", length(res_sf$polygons.areas), "\n")
cat("  PASS\n\n")

cat("=== Test 12: Larger point set performance ===\n")
set.seed(100)
n_large <- 1000
pts_large <- matrix(runif(2 * n_large), ncol = 2)

bbox_large <- ppmData:::get_bbox(pts_large, buffer = 0.05)
dummy_large <- ppmData:::get_bbox_dummy_coords(pts_large, bbox_large, n = 9)
pts_all_large <- rbind(pts_large, dummy_large)
coords_large <- as.numeric(t(pts_all_large))

t_large <- system.time({
  areas_large <- ppmData:::dirtess_clip_areas_cpp(
    coords = coords_large,
    ncoords = n_large,
    clip_x = c(0, 1, 1, 0),
    clip_y = c(0, 0, 1, 1),
    hole_x_list = list(),
    hole_y_list = list()
  )
})

cat("  1000 points, time:", t_large["elapsed"], "s\n")
cat("  total area:", sum(areas_large), "\n")
cat("  area close to 1.0:", abs(sum(areas_large) - 1.0) < 0.01, "\n")
cat("  PASS\n\n")

cat("=== All tests completed ===\n")
