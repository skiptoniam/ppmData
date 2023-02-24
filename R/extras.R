expectClass <- function (x, expected_classes) {

  cl <- class(x)
  if (!inherits(x,expected_classes)) {
    text <- sprintf('This object should be of class: %s, but had class%s: %s',
                    paste(expected_classes, collapse = ', '),
                    ifelse(length(cl) > 1, 'es', ''),
                    paste(cl, collapse = ', '))

    stop(text)
  }
}

list2numeric <- function(x){
  as.numeric(as.matrix(x))
}


