(in-package regression)

(defun float-sigma (x)
  "Calculate logistic function x -> σ(x)

        1
σ = ---------
    1+exp(-x)
"
  (declare (single-float x)
	   (optimize speed))
  (/ 1s0 (+ 1s0 (exp (- x)))))

(defun float-dsigma (diff sigma)
  (* sigma (- 1s0 sigma) diff))

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

(defun make-random-array (x y scale)
  "Make random array of single floats in scale."
  (let ((res
	  (make-array (list x y) :element-type 'single-float)))
    (dotimes (i  x res)
      (dotimes (j y)
	(setf (aref res i j) (- scale (random (* 2 scale))))))))
