(defpackage "LINEAR-ALGEBRA"
  (:use #:cl #:rt)
  (:export "TIMES" "TIMES-INTO"
	   "LINEAR-COMBINATION" "LINEAR-COMBINATION-INTO" "LINEAR-UPDATE"
	   "TIMES-TRANSPOSED" "TIMES-TRANSPOSED-INTO"
	   "TRACE-TIMES-TRANSPOSED"
	   "TIMES-REV-TRANSPOSED" "M-"))

(defpackage #:regression
   (:use #:cl #:rt #:linear-algebra)
   (:export "LINEAR-PROPAGATE" "LOGISTIC-PROPAGATE"
	    "LINEAR-REGRESSION-ITERATION" "LOGISTIC-REGRESSION-ITERATION"
	    "MAKE-RANDOM-ARRAY")
   (:documentation "Regression tools."))
