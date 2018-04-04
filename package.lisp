(defpackage "LINEAR-ALGEBRA"
  (:use #:cl #:rt)
  (:export "TIMES" "TIMES-INTO"
	   "LINEAR-COMBINATION"
	   "TIMES-TRANSPOSED" "TIMES-TRANSPOSED-INTO"))

(defpackage #:regression
   (:use #:cl #:rt #:linear-algebra)
   (:export "LINEAR-PROPAGATE" "LOGISTIC-PROPAGATE"
	    "LINEAR-REGRESSION-ITERATION" "LOGISTIC-REGRESSION-ITERATION"
	    "MAKE-RANDOM-ARRAY")
   (:documentation "Regression tools."))
