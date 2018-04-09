(in-package regression)

;; neural network: 10 input cells, 3 middle cells, 5 output cells
(defvar *a* (make-random-array 10 3 1.0))
(defvar *b* (make-random-array 3 5 1.0))

(defvar *guess-a* (make-random-array 10 3 1.0))
(defvar *guess-b* (make-random-array 3 5 1.0))


(defun true-predict (x)
  (times (apply-fn #'float-sigma
		   (times x *a*)) 
	 *b*))

(defvar *xses* (make-random-array 50 10 1s0))
(defvar *yses* (true-predict *xses*))

(defun forward (A B)
  (let* ((middle (apply-fn #'float-sigma (times *xses* A)))
	 (final (times middle B)))
    (values middle final)))

(defun backward (middle final y A B &optional (sigma 0.3) (rho 1s0))
  ;; copy from logistic-regression-iteration
  (let* ((my-y final)
	 (dy-second-layer
	   (linear-combination 1s0
			       (copy-array my-y)
			       -1s0 y))
	 (f (times-transposed dy-second-layer dy-second-layer))

	 (b-diff (times  dy-second-layer middle))
	 (mid-diff (times B dy-second-layer))
	 (a-diff (apply-fn2 #'float-dsigma
			    middle mid-diff)))
    (linear-combination rho A sigma a-diff)
    (linear-combination rho B sigma b-diff)
    (values f A B)))
