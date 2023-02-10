BEGIN {
  print "["
  first_line = 1
}

{
  if (first_line == 1) {
    print $0
    first_line = 0
  } else {
    print "," $0
  }
}

END {
  print "]"
}