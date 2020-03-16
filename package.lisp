(cz.zellerin.doc:defpackage "LINEAR-ALGEBRA"
  (:use #:cl #:rt)
  (:export "WITH-MATRIXES" "TRANSPOSE"
	   "TRACE" "SCALAR"

	   "*MATRIX-FIELD*" "MATRIX-OPTIMIZE")
  (:sections @matrix-ops
	     @variables
	     @matrix-ops-implementation)
  (:documentation
   "For the regression, I need linear combinations and matrix (including
row and column vector) multiplication. I met quite a few LA libraries,
some native, some using FFI, and most of them is a bit too
complicated. Also, I wanted to get a feeling on speed of simple code."))

(cz.zellerin.doc:defpackage #:regression
   (:use #:cl #:rt #:linear-algebra)
   (:export "LINEAR-PROPAGATE" "LOGISTIC-PROPAGATE"
	    "LINEAR-REGRESSION-ITERATION" "LOGISTIC-REGRESSION-ITERATION"
	    "MAKE-RANDOM-ARRAY")
   (:documentation "Regression tools.")
   (:sections @utilities))
