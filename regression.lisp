(in-package regression)

;  "APPLY-FN" "APPLY-FN2"

;;; * Linear regression
;;; model: y'=xA
;;; Square error: F = Tr (y-y')(y-y')†
;;; function: F+½ρ² Tr AA†
;;; Variation against A: 2(y'-y)†x + ρ²A

;; test model
(defparameter *true-a* (make-random-array 3 1 1s0))
(defparameter *x* (make-random-array 10 3 5s0))
(defparameter *y* (times *more-x* *true-a*))
(defparameter *logistic-y* (apply-fn #'float-sigma (copy-array *y*)))
(defparameter *test-A* (make-random-array 3 1 1s2))

(defun linear-regression-iteration (y A x sigma rho)
  "Iteration step; modifies A."
  (let ((dy (linear-combination 1s0 (times x A)
				-1s0 y)))
    (values
     (times-transposed dy dy)
     dy
     (linear-combination rho A sigma (times-transposed x dy ))
     *true-a*)))

(defun logistic-regression-iteration (y A x sigma rho)
  "Iteration step; modifies A to get closer."
  (let ((dy (linear-combination 1s0
				(apply-fn #'float-sigma
					  (times x A))
				-1s0 y)))
    (values
     (times-transposed dy dy)
     dy
     (linear-combination rho A sigma
			 (times-transposed x (apply-fn2 #'float-dsigma dy y)))
     *true-a*)))

(defun check-linear (count)
  (let ((*test-A* (make-random-array 3 1 1s0)))
    (dotimes (i count)
      (print (linear-regression-iteration *y* *test-A* *more-x* -6s-3 1s0)))))

(defun check-logistic (count)
  (let ((*test-A* (make-random-array 3 1 1s0)))
    (dotimes (i count)
      (print (logistic-regression-iteration *logistic-y* *test-A* *more-x* 0.3 1s0)))))

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

; (try-sigmas-logistic *logistic-y* *test-A* *more-x* 1s0 1s0)

(defun try-sigmas-linear (&rest pars)
  "Try range of sigmas (two orders) around base for linear regression."
  (apply 'try-sigmas 'linear-regression-iteration pars))

; (try-sigmas-linear *y* *test-A* *more-x* -1s-2 1s0)
