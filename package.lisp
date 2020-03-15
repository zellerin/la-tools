(cz.zellerin.doc:defpackage "LINEAR-ALGEBRA"
  (:use #:cl #:rt)
  (:export "WITH-MATRIXES" "TRANSPOSE"
	   "TRACE" "SCALAR"

	   "*MATRIX-FIELD*" "MATRIX-OPTIMIZE")
  (:sections @matrix-ops))

(cz.zellerin.doc:defpackage #:regression
   (:use #:cl #:rt #:linear-algebra)
   (:export "LINEAR-PROPAGATE" "LOGISTIC-PROPAGATE"
	    "LINEAR-REGRESSION-ITERATION" "LOGISTIC-REGRESSION-ITERATION"
	    "MAKE-RANDOM-ARRAY")
   (:documentation "Regression tools.")
   (:sections @utilities))
