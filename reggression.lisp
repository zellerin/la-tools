(defpackage #:regression
   (:use #:cl #:linear-algebra)
   (:export "LINEAR-PROPAGATE" "LOGISTIC-PROPAGATE"
	    "LINEAR-REGRESSION-ITERATION" "LOGISTIC-REGRESSION-ITERATION"
	    "MAKE-RABDOM-ARRAY")
   (:documentation "lorem ipsum"))

(in-package regression)

;  "APPLY-FN" "APPLY-FN2"

;;; * Linear regression
;;; model: y'=xA
;;; Square error: F = Tr (y-y')(y-y')†
;;; function: F+½ρ² Tr AA†
;;; Variation against A: 2(y'-y)†x + ρ²A

(defun make-random-array (x y scale)
  (let ((res
	  (make-array (list x y) :element-type 'single-float)))
    (dotimes (i  x res)
      (dotimes (j y)
	(setf (aref res i j) (random scale))))))

;; test model
(defparameter *true-a* (make-random-array 3 1 1s0))

(defparameter *more-x* (make-random-array 10 3 5s0))

(defparameter *y* (times *more-x* *true-a*))

(defparameter *test-A* (make-random-array 3 1 1s0))

(defun linear-regression-iteration (y A x sigma rho)
  "Iteration step; modifies A."
  (let ((dy (linear-combination 1s0 (times x A)
				-1s0 y)))
    (values
     (times-transposed dy dy)
     dy
     (linear-combination rho A sigma (times-transposed x dy ))
     *true-a*)))

(defun apply-fn (fn A)
  (let ((rows (array-dimension A 0))
	(cols (array-dimension A 1)))
  (declare
   ((simple-array single-float) A)
   (optimize speed)
   (compiled-function fn))
    (dotimes (i rows A)
      (dotimes (j cols)
	(setf (aref A i j) (funcall fn (aref A i j)))))))

(defun apply-fn2 (fn A B)
  (let ((rows (array-dimension A 0))
	(cols (array-dimension A 1)))
  (declare
   ((simple-array single-float) A B)
   (optimize speed)
   (compiled-function fn))
    (dotimes (i rows A)
      (dotimes (j cols)
	(setf (aref A i j) (funcall fn (aref A i j) (aref B i j)))))))

(defun float-sigma (x)
  (declare (single-float x)
	   (optimize speed))
  (/ 1s0 (+ 1s0 (exp x))))

(defparameter *logistic-y* (apply-fn #'float-sigma (copy-array *y*)))

(defun float-dsigma (diff sigma)
  (* sigma (- 1s0 sigma) diff))

(defun logistic-regression-iteration (y A x sigma rho)
  "Iteration step; modifies A."
  (let ((dy (linear-combination 1s0
				(apply-fn #'float-sigma
					  (times x A))
				-1s0 y)))
    (values
     (times-transposed dy dy)
     dy
     (linear-combination rho A sigma
			 (times-transposed x (apply-fn2 #'float-dsigma dy y)))
     *true-a*
    )))


(defun check-linear (count)
  (let ((*test-A* (make-random-array 3 1 1s0)))
    (dotimes (i count)
      (print (linear-regression-iteration *y* *test-A* *more-x* -5s-3 1s0)))))

(defun check-logistic (count)
  (let ((*test-A* (make-random-array 3 1 1s0)))
    (dotimes (i count)
      (print (logistic-regression-iteration *logistic-y* *test-A* *more-x* 0.3 1s0)))))

(defun copy-array (a)
  (let* ((type (array-element-type a))
	 (res
	   (make-array (array-dimensions a)
		       :element-type type)))
    (map-into (make-array (apply #'* (array-dimensions a))
	       :element-type type
	       :displaced-to res)
	      'identity
	      (make-array (apply #'* (array-dimensions a))
			  :element-type type
			  :displaced-to a))
    res))

(defun try-sigmas (fn y fixed-A x base rho)
  (loop for i from -1s0 to 1s0 by 0.1
	for sigma = (* base (expt 10s0 i))
	do
	   (let ((A (copy-array fixed-A)))
	     (funcall fn y a x sigma rho)
		       (print (list sigma (aref
					   (funcall fn y a x sigma rho)
					   0 0))))))

(defun try-sigmas-logistic (&rest pars)
  (apply 'try-sigmas 'logistic-regression-iteration pars))

; (try-sigmas-logistic *logistic-y* *test-A* *more-x* 1s0 1s0)

(defun try-sigmas-linear (&rest pars)
  (apply 'try-sigmas 'linear-regression-iteration pars))

; (try-sigmas-linear *y* *test-A* *more-x* -1s-2 1s0)

