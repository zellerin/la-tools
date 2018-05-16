;;; -*- mode: lisp -*-

;;; Note: this is regression in the mathematical sense (linear
;;; regression and logistic regression - fitting parameters to a
;;; model), not in the software sense (regression tests)
(in-package regression)

(defun M- (a b)
  (linear-combination 1s0 a -1s0 b))

(declaim (inline logistic-estimate linear-estimate
		 logistic-grad-A linear-grad-A))

(defun logistic-estimate (X A)
  "Logistic regression estimate given coefficients A and independent variables X."
  (apply-fn #'float-sigma (times X A)))

(defun linear-estimate (X A)
  "Logistic regression estimate given coefficients A and independent variables X."
  (times X A))

(defun logistic-grad-A (X err Y)
  (times-transposed X (apply-fn2 #'float-dsigma err Y)))

(defun linear-grad-A (X err Y)
  (declare (ignore Y))
  (times-transposed X err))

(eval-when (:compile-toplevel)
  (defvar *pairs* '(estimate grad-a)))

(defmacro define-pair (name pars
		       docstring-fmt
		       &body body)
  `(progn
     (eval-when (:compile-toplevel)
       (pushnew ',name *pairs*))
     (defun ,name ,pars
       ,docstring-fmt
       (declare (ignorable ,@(mapcar
			      ;; cannot use ignore, as parameter may
			      ;; be used in default calculation for
			      ;; optional/key parameter
			      (lambda (a) (if (consp a) (car a) a))
			   (set-difference pars '(&aux &rest &key &optional)))))
       (error "This function cannot be called outside DEFINE-PAIR"))
     ,@(loop for prefix in '("LINEAR" "LOGISTIC")
	     collect
	     (labels ((@ (name)
			(alexandria:symbolicate prefix "-" name))
		      (flet-list-item (name)
			`(,name (&rest pars) (apply #',(@ name) pars))))

	       `(defun ,(@ name) ,pars
		  ,(format nil docstring-fmt prefix)
		  (flet ,(mapcar #'flet-list-item *pairs*)
		    ,@body))))))

(define-pair regression-iteration (Y A X sigma rho)
  "Update matrix with ~(~A~) regression coeficients to
better match observed data.

Takes as parameter
- the model data original regression coefficients matrix A, matrix of
  the independent values X, matrix of observed values Y
- the regression process data: rho below 1s0 prevents overfitting, and
  sigma determines the speed of the regression.

Return three values: square of the error vector (difference between
estimated and provided Y), error vector itself and updated matrix with
the regression coefficients."
  (let* ((estimated-Y (estimate X A))
	 (err (M- estimated-Y y))
	 (a-diff (grad-A X (copy-array err) estimated-Y)))
    (values
	    (trace-times-transposed err err)
	    err (linear-update rho A sigma a-diff))))

(define-pair regression-iterations (Y A X sigma alpha count
				      &optional out (sampling 20)
				      &aux (rho (+ 1s0 (* alpha sigma))))
    "Run COUNT iterations, optionally logging error(s)."
  (dotimes (i count)
    (let ((err (regression-iteration y a x sigma rho)))
      (when (and out (zerop (mod i sampling)))
	(let ((A-err (* alpha 2s0 (trace-times-transposed A A))))
	  (format out "~a ~a ~a~%"
		  err A-err (+ err A-err)))))))

(define-pair test-case (samples indeps &optional (deps 1))
    "Make a test case for ~(~A~) regression.

Returns X, Y and initial A as values."
  (let ((X (make-random-array samples indeps 1s0)))
    (values X (estimate X (make-random-array indeps deps 1s0))
	    (make-random-array indeps deps 1s-5))))

(defmacro with-test-case ((samples indeps deps) &body body)
  `(multiple-value-bind (X Y A) (test-case ,samples ,indeps ,deps)
    ,@body))

(define-pair try-sigmas (y raw-x base &key
			   (fixed-A (make-random-array (array-dimension raw-x 1) 1 2s-2))
			   (count 1000)
			   (alpha 0s0)
			   (sampling 30)
			   out)
  "Try range of sigmas (two orders) around base."
  (multiple-value-bind (x norm) (normalize raw-x)
    (declare (ignore norm))
    (loop for i from -1s0 to 1s0 by 0.5
	  for sigma = (* base (expt 10s0 i))
	  for rho = (- 1s0 (* alpha sigma))
	  do (format out "~2%\"~a\"~%" sigma)
	  collect
	     (let ((A (copy-array fixed-A)))
	       (regression-iterations y a x sigma alpha count out sampling)))))

(define-pair check-regression (count &optional
				     (sampling 20)
				     (sigma -6s-3))
  "Test COUNT rounds of regression on random matrixes"
  (with-test-case (100 120 1)
    (try-sigmas Y X sigma)
    (dotimes (i count)
      (let ((F (regression-iterations Y A X
				      sigma 0s0
				      out sampling)))))))

(define-pair get-coefficients (y raw-x &key
				 (A (make-random-array (array-dimension raw-x 1) 1 2s-2))
				 (sigma (- (/ 1s0 (array-dimension y 0)))) (alpha 0s0)
				 (rho (+ 1s0 (* alpha sigma))))
    "Get coefficients for ~(~A~) regression from data.

Assumes that first row of X is all ones. This corresponds to getting
linear term, and it is also used for normalization purposes."
  (multiple-value-bind (x norm) (normalize raw-x)
    (dotimes (i 4000)
      (regression-iteration y A x sigma rho))
    (make-array (array-dimension A 0)
		:element-type 'single-float
		:displaced-to (times norm A))))
