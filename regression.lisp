(in-package regression)

;;; * Linear regression
;;; model: y'=xA
;;; Square error: F = Tr (y-y')(y-y')†
;;; function: F+½ρ² Tr AA†
;;; Variation against A: 2(y'-y)†x + ρ²A

(in-package regression)

(defun M- (a b)
  (linear-combination 1s0 a -1s0 b))

(defun times-sigma (X A)
  "Logistic regression estimate given coefficients A and independent variables X."
  (apply-fn #'float-sigma (times X A)))

(defun logistic-grad-A (X err Y)
   (times-transposed X (apply-fn2 #'float-dsigma err Y)))

(defun l-regression-iteration (Y A X sigma rho estimate-fn grad-A-fn)
  "Update matrix with linear or logistic regression coeficients to
better match observed data.

Takes as parameter 
- the model data original regression coefficients matrix A, matrix of
  the independent values X, matrix of observed values Y
- the regression process data: rho below 1s0 prevents overfitting, and
  sigma determines the speed of the regression.
- the regression model data (functions specific for linear or logistic regression)

Return three values: square of the error vector (difference between
estimated and provided Y), error vector itself and updated matrix with
the regression coefficients."
  (let* ((estimated-Y (funcall estimate-fn X A))
	 (err (M- estimated-Y y))
	 (a-diff (funcall grad-A-fn X err estimated-Y)))
    (values (times-transposed err err)
	    err (linear-update rho A sigma a-diff))))

(defun logistic-regression-iteration (Y A X sigma rho)
  "See l-regression-iteration for details."
  (l-regression-iteration Y A X sigma rho #'times-sigma #'logistic-grad-A))

(defun linear-regression-iteration (Y A X sigma rho)
  "See l-regression-iteration for details."
  (l-regression-iteration Y A X sigma rho #'times
			  (lambda (X err A)
			    (declare (ignore A))
			    (times-transposed X err))))

(defun check-linear (count &optional (sampling 20))
  (let ((*test-A* (make-random-array 3 1 1s0)))
    (dotimes (i count)
      (let ((F
	      (linear-regression-iteration *y* *test-A* *x* -6s-3 1s0)))
	(when (zerop (mod i sampling))
	  (print (aref F 0 0)))))))

(defun check-logistic (count pars)
  (let ((*test-A* (make-random-array 3 1 1s0)))
    (dotimes (i count)
      (print (apply 'logistic-regression-iteration pars)))))

(defun try-sigmas (fn y fixed-A x base rho &optional (count 10))
  "Try range of sigmas (two orders) around base."
  (loop for i from -1s0 to 1s0 by 0.1
	for sigma = (* base (expt 10s0 i))
	do
	   (let ((A (copy-array fixed-A)))
	     (princ sigma)
	     (dotimes (i count)
	       (format t " ~10e " (aref (funcall fn y a x sigma rho) 0 0)))
	     (terpri))))

(defun try-sigmas-logistic (&rest pars)
  "Try range of sigmas (two orders) around base for logistic regression."
  (apply 'try-sigmas 'logistic-regression-iteration pars))

; (try-sigmas-logistic *logistic-y* *test-A* *x* -0.006 1s0)

(defun try-sigmas-linear (&rest pars)
  "Try range of sigmas (two orders) around base for linear regression."
  (apply 'try-sigmas 'linear-regression-iteration pars))

; (try-sigmas-linear *y* *test-A* *x* -1s-2 1s0)
