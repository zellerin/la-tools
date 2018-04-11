(in-package regression)

;; neural network: 10 input cells, 3 middle cells, 5 output cells
(defvar *a* (make-random-array 10 3 1.0))
(defvar *b* (make-random-array 3 1 1.0)) ; hidden final

(defvar *guess-a* (make-random-array 10 3 1.0))
(defvar *guess-b* (make-random-array 3 1 1.0))


(defun true-predict (x)
  (times (apply-fn #'float-sigma
		   (times x *a*))
	 *b*))

(defvar *xses* (make-random-array 50 10 1s0))
(defvar *yses* (true-predict *xses*))

(defun forward-layer (x A)
  (apply-fn #'float-sigma (times x A)) ; cases hidden
)

(defun backward-linear-layer (dF/dy x A)
  ;; y        cases  y-size
  ;; dF/dy    cases y-sizes
  ;; A        x-size y-size
  ;; x        cases  x-size

  ;; dF/dA,   x-size y-size
  ;; dF/dx+   x-size cases
  (values (times-transposed x dF/dy)
	  (times-rev-transposed dF/dy A)))

(defun backwars-logistic-layer (y dF/dy x A)
  (backward-linear-layer (apply-fn2 #'float-dsigma dF/dy y)
			 x A))

(defun forward (A B)
  (let* ((middle (apply-fn #'float-sigma (times *xses* A))) ; cases hidden
	 (final (times middle B))) ; cases 1
    (values middle final)))

(defun backward (x hidden y A B &optional (a-sigma -0.3) (b-sigma -0.3) (rho 1s0))
  ;; copy from logistic-regression-iteration
  (let ((dy (linear-combination -1s0 y 1s0 *yses*)))
    (multiple-value-bind (var-B dF/dhidden) (backward-linear-layer dy hidden B)
      (multiple-value-bind (var-A dF/dx) (backwars-logistic-layer hidden dF/dhidden x A)
	(declare (ignore dF/dx))
	(linear-update rho A a-sigma var-A)
	(linear-update rho B b-sigma var-B)
	(times-transposed dy dy)))))
