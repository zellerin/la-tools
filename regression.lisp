(in-package regression)

(defmacro do-regression (estimator updater
			 (Y A sigma rho count)
			 &body body)
  "Update up to COUNT times matrix A, and execute BODY in each iteration.

In each step, A is set to ρA+σU, where U is provided updater.

Parameter ρ below 1s0 prevents overfitting, and σ determines the speed of the regression.

The body is executed with these bindings:
- ESTIMATE for estimated value,
- ERR for difference between the ESTIMATE and Y
- Iteration count I
"
  `(loop
     with err = (make-array (array-dimensions Y)
			    :element-type 'single-float)
     and estimate = (make-array (array-dimensions Y)
				:element-type 'single-float)
     for i from 0 to ,count
     do
	(with-matrixes ,estimator :declarations nil
	  :target estimate)
	(with-matrixes (- estimate ,Y)
	  :declarations nil
	  :target err)
	(with-matrixes (+ (* ,rho ,A)
			  (* ,sigma ,updater))
	  :target ,A
	  :declarations ((scalar ,rho ,sigma)))
	(setf ,sigma (* ,sigma 1.0001))
	(progn ,@body)))

(defun linear-regression-iterations (Y A X msigma alpha count out
				     &aux
				       (sigma (/ msigma (array-dimension X 0)))
				       (rho (+ 1s0 (* alpha msigma))))
  (declare (single-float sigma rho)
	   ((simple-array single-float) Y A X))
  (do-regression (* X A) (* (transpose X) err)
      (Y A sigma rho count)
    (when (and out (zerop (mod i 100)))
      (print (trace-times-transposed err err) out))))

(defun logistic-regression-iterations (Y A X msigma alpha count out
				     &aux
				       (sigma (/ msigma (array-dimension X 0)))
				       (rho (+ 1s0 (* alpha msigma))))
  (declare (single-float sigma rho)
	   ((simple-array single-float) Y A X))
  (let ((samples (array-dimension X 0))
	(parameters (array-dimension X 1))
	(dependents (array-dimension A 1)))
    (when out (format out "~&# Samples: ~d Parameters: ~d Dependends: ~d~%" samples parameters dependents))
    (do-regression (map sigma (* X A))
	(* (transpose X) (map dsigma err estimate))
	(Y A sigma rho count)
      (when (and out (zerop (mod i 100)))
	(let ((e (/ (trace-times-transposed err err) 2.0 samples dependents))
	      (a-err (/ (* (- 1 rho) (trace-times-transposed A A)) 2.0 dependents)))
	  (format out "~a ~a ~a~%"
		  e A-err  (+ e A-err)))))))

(defun get-coefficients (fn y raw-x &key
				   (A (make-random-array (array-dimension raw-x 1) 1 2s-2))
				   (sigma (- (/ 1s0 (array-dimension y 0))))
				   (alpha 0s0)
				   out)
    "Get coefficients for linear regression from data.

Assumes that first row of X is all ones. This corresponds to getting
linear term, and it is also used for normalization purposes."
  (multiple-value-bind (x norm) (normalize raw-x)
    (funcall fn y A x sigma alpha 4000 out)
    (make-array (array-dimension A 0)
		:element-type 'single-float
		:displaced-to (with-matrixes (* norm A)))))

(defun linear-get-coefficients (&rest args)
  (apply #'get-coefficients #'linear-regression-iterations args))

(defun logistic-get-coefficients (&rest args)
  (apply #'get-coefficients #'logistic-regression-iterations args))

(defmacro with-test-case ((samples indeps deps) &body body)
  `(multiple-value-bind (X Y A) (test-case ,samples ,indeps ,deps)
    ,@body))

(defun linear-test-case (samples indeps &optional (deps 1))
    "Make a test case for linear regression.

Returns X, Y and initial A as values."
  (let ((X (make-random-array samples indeps 1s0)))
    (values X (linear-estimate X (make-random-array indeps deps 1s0))
	    (make-random-array indeps deps 1s-5))))

#+nil (defun linear-check-regression (count &key
				     (samples 10)
				     (indeps 10)
				     (deps 1)
				     (sampling 20)
				     (sigma -6s-3)
				     (alpha 0s0)
				     (out *standard-output*))
  "Test COUNT rounds of regression on random matrixes"
  (multiple-value-bind (X Y A) (linear-test-case samples indeps deps)
    (linear-regression-iterations Y A X
				  sigma alpha
				  count out sampling)))
