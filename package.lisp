(defpackage "LINEAR-ALGEBRA"
  (:use #:cl #:rt)
  (:export "TIMES" "TIMES-INTO"
	   "LINEAR-COMBINATION" "LINEAR-COMBINATION-INTO" "LINEAR-UPDATE"
	   "TIMES-TRANSPOSED" "TIMES-TRANSPOSED-INTO"
	   "TIMES-REV-TRANSPOSED"))

(defpackage #:regression
   (:use #:cl #:rt #:linear-algebra)
   (:export "LINEAR-PROPAGATE" "LOGISTIC-PROPAGATE"
	    "LINEAR-REGRESSION-ITERATION" "LOGISTIC-REGRESSION-ITERATION"
	    "MAKE-RANDOM-ARRAY")
   (:documentation "Regression tools."))
