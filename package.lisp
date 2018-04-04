(defpackage #:regression
   (:use #:cl)
   (:export "LINEAR-PROPAGATE" "LOGISTIC-PROPAGATE"
	    "LINEAR-REGRESSION-ITERATION" "LOGISTIC-REGRESSION-ITERATION"
	    "MAKE-RANDOM-ARRAY")
   (:documentation "Regression tools."))
