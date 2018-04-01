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
(defvar *true-a* (make-random-array 3 1 1s0))

(defvar *more-x* (make-random-array 10 3 5s0))

(defvar *y* (times *more-x* *true-a*))

(defvar *test-A* (make-random-array 3 1 1s0))

(defvar *dy* (linear-combination 1s0 (times *more-x* *test-A*)
				 -1s0 *y*))

(defun linear-regression-iteration (y A x sigma rho)
  "Iteration step; modifies A."
  (let ((dy (linear-combination 1s0 (times x A)
				-1s0 y)))
    (values
     dy
     (linear-combination rho A sigma (times-transposed x dy ))
     *true-a*
     (times-transposed dy dy))))

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

(defvar *logistic-y* (apply-fn #'float-sigma *y*))

(defun float-dsigma (diff sigma)
  (* sigma (- 1s0 sigma) diff))

(defun logistic-regression-iteration (y A x sigma rho)
  "Iteration step; modifies A."
  (let ((dy (linear-combination 1s0
				(apply-fn #'float-sigma
					  (times x A))
				-1s0 y)))
    (values
     dy
     (linear-combination rho A sigma
			 (times-transposed x (apply-fn2 #'float-dsigma dy y)))
     *true-a*
     (times-transposed dy dy))))


(defun check-linear (count)
  (let ((*test-A* (make-random-array 3 1 1s0)))
    (dotimes (i count (linear-regression-iteration *y* *test-A* *more-x* -1s-4 1s0))
      (linear-regression-iteration *y* *test-A* *more-x* -1s-4 1s0))))

(defun copy-array (a)
  (make-array (array-dimensions a)
	      :element-type (array-element-type a)
	      :initial-contents a))
