(in-package regression)

(macrolet ((make-updater-step (name estimator updater)
	     `(defun ,name (A Y X estimate err sigma rho)
		(declare ((simple-array single-float) Y A X estimate err))
		(declare (single-float rho sigma))
		(with-matrixes ,estimator :declarations nil :target estimate)
		(with-matrixes (- Y estimate) :declarations nil :target err)
		(with-matrixes (+ (* rho A) (* sigma ,updater)) :target A
		  :declarations ((scalar rho sigma))))))

  (make-updater-step linear-updater (* X A) (* (transpose X) err))
  (make-updater-step logistic-updater
		     (map sigma (* X A))
		     (* (transpose X) (map dsigma err estimate))))

(macrolet ((make-updater-step (name estimator updater)
	     `(defun ,name (A Y X estimate err sigma rho)
		(declare ((simple-array single-float) Y A X estimate err))
		(declare (single-float rho sigma))
		(with-matrixes ,estimator :declarations nil :target estimate)
		(with-matrixes (- Y estimate) :declarations nil :target err)
		(with-matrixes (+ (map (lambda (a)
					 (cond ((> a rho) (- a rho))
					       ((> (- rho) a) (+ rho a))
					       (t 0))) A)
				  (* sigma ,updater)) :target A
		  :declarations ((scalar rho sigma))))))

  (make-updater-step linear-updater-l1reg (* X A) (* (transpose X) err))
  (make-updater-step logistic-updater-l1reg
		     (map sigma (* X A))
		     (* (transpose X) (map dsigma err estimate))))

(defun do-regression-fn (update-fn
			 Y A X sigma rho count
			 body-fn)
  "Update up to COUNT times matrix A, and execute BODY in each iteration.

In each step, A is set to ρA+σU, where U is provided updater.

Parameter ρ below 1s0 prevents overfitting, and σ determines the speed of the regression.
They are relate

The ESTIMATOR is a matrix expression passed to `with-matrixes' macro.

The body is executed with these bindings:
- ESTIMATE for estimated value,
- ERR for difference between the ESTIMATE and Y
- Iteration count I
and can use (RETURN) to terminate the loop, or change the parameters.

ESTIMATOR and ERR are calculated even when not needed. This may change in future.
"
  (loop
    with err = (make-array (array-dimensions Y)
			   :element-type 'single-float)
    and estimate = (make-array (array-dimensions Y)
			       :element-type 'single-float)
    for i from 0 to count
    do
       (funcall update-fn A Y X estimate err sigma rho)
       (funcall body-fn i err)))

(defun do-normalized-iteration (update-fn
				   Y A X iteration-speed regularization-cost-per-item count out
				   tracing)
  (let* ((samples (array-dimension X 0))
	 (parameters (array-dimension X 1))
	 (dependents (array-dimension A 1))
	 (sigma (/ iteration-speed samples dependents))
	 (rho (- 1s0
		  (/
		   (* regularization-cost-per-item iteration-speed)
		   parameters dependents))))
    (assert (>= 1.0 rho 0.0))
    (assert (>= sigma 0.0))
    (when out
      (format out "~&# Samples: ~d Parameters: ~d Dependends: ~d~%" samples parameters dependents)
      (format out "~&# rho: ~f sigma: ~s~%" rho sigma))
    (do-regression-fn update-fn
      Y A X sigma rho count
      (lambda (i err)
	(when (and out (zerop (mod i tracing)))
	  (let ((e (/ (trace-times-transposed err err) 2.0 samples dependents))
		(a-err (/ (* regularization-cost-per-item
			     (trace-times-transposed A A))
			  dependents parameters 2.0)))
	    (format out "~a ~a ~a ~a~%"
		    i e A-err (+ e A-err))))))))

(defun linear-regression-iterations (Y A X msigma alpha count out tracing)
  (declare ((simple-array single-float) Y A X))
  (do-normalized-iteration
      #'linear-updater Y A X msigma alpha count out tracing))

(defun logistic-regression-iterations (Y A X msigma alpha count out tracing)
  (do-normalized-iteration #'logistic-updater
    Y A X msigma alpha count out  tracing))

(defun get-coefficients (fn y raw-x &key
				      (A (make-random-array (array-dimension raw-x 1) 1 1s-2))
				      (sigma (- (/ 1s0 (array-dimension y 0))))
				      (alpha 0s0)
				      out (tracing 100)
				      (count 1000))
  "Get coefficients for linear regression from data.

Assumes that first row of X is all ones. This corresponds to getting
linear term, and it is also used for normalization purposes."
  (multiple-value-bind (x norm) (normalize raw-x)
    (do-normalized-iteration fn y A x sigma alpha count out tracing)
    (make-array (array-dimension A 0)
		:element-type 'single-float
		:displaced-to (with-matrixes (* norm A)
				:optimize ((speed 1))))))

(defun linear-get-coefficients (&rest args)
  (apply #'get-coefficients #'linear-regression-iterations args))

(defun logistic-get-coefficients (&rest args)
  (apply #'get-coefficients #'logistic-regression-iterations args))

(defmacro with-test-case ((samples indeps deps) &body body)
  `(multiple-value-bind (X Y A) (test-case ,samples ,indeps ,deps)
     ,@body))
