;;; model: y'=xA
;;; Square error: F = Tr (y-y')(y-y')†
;;; function: F+½ρ² Tr AA†
;;; Variation against A: 2(y'-y)†x + ρ²A

(defun random-array (x y scale)
  (let ((res
	  (make-array (list x y) :element-type 'single-float)))
    (dotimes (i  x res)
      (dotimes (j y)
	(setf (aref res i j) (random scale))))))

;; test model
(defvar *true-a* (random-array 3 1 1s0))

(defvar *more-x* (random-array 10 3 5s0))

(defvar *y* (times *more-x* *true-a*))

(defvar *test-A* (random-array 3 1 1s0))

(defvar *dy* (linear-combination 1s0 (times *more-x* *test-A*)
				 -1s0 *y*))

(defvar *dA* (times-transposed *more-x* *dy*))

(defun advance (y A x sigma rho)
  (let ((dy (linear-combination 1s0 (times x A)
				-1s0 y)))
    (values
     dy
     (linear-combination rho A sigma (times-transposed x dy))
     *true-a*
     (times-transposed dy dy))))
