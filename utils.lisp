(in-package regression)

(defun float-sigma (x)
  "Calculate logistic function x -> σ(x)

        1
σ = ---------
    1+exp(-x)

for a single-float parameter X. Care is taken not to underflow for
large negative values of x."
  (declare (single-float x)
	   (optimize speed))
  (if (> x #.(- 1s0 (log most-positive-single-float)))
      (/ 1s0 (+ 1s0 (exp (- x))))
      0s0))

(defun float-dsigma (diff sigma)
  "Help calculate derivation of function of sigma, using formula

f(sigma(x))' = f'(sigma) . sigma . (1-sigma)."
  (* sigma (- 1s0 sigma) diff))

(defun apply-fn (fn A)
  "Destructively apply function fn to all cells of a 2D matrix A."
  (let ((rows (array-dimension A 0))
	(cols (array-dimension A 1)))
  (declare
   ((simple-array single-float) A)
   (optimize speed)
   (compiled-function fn))
    (dotimes (i rows A)
      (dotimes (j cols)
	(setf (aref A i j) (funcall fn (aref A i j)))))))

(defun apply-fn2 (fn A B )
  "Return result of application of function fn to all respective cells of a 2D
matrix A and B (same dimensions assumed) and store result to cells in A.."
  (let* ((rows (array-dimension A 0))
	 (cols (array-dimension A 1))
	 (res (make-array (list rows cols)
			  :element-type (array-element-type A))))
  (declare
   ((simple-array single-float) A B res)
   (optimize speed)
   (compiled-function fn))
    (dotimes (i rows A)
      (dotimes (j cols)
	(setf (aref res i j) (funcall fn (aref A i j) (aref B i j)))))))


(defun copy-array (a)
  "Create a copy of an array with same dimensions and element type."
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

(defun normalize (x)
  "Assuming that first column of X is all ones, return normalized matrix X and normalizer A, i.e.,
X -> (X' A) where
- X' = X.A,
- and average of each X column is 0 and standard
deviation 1."
  (let* ((correlation (times-transposed X X))
	 (size (array-dimension X 1))
	 (A (make-array (list size size) :element-type 'single-float
					 :initial-element 0s0)))
    (setf (aref A 0 0) 1s0)
    (loop with len = (array-dimension X 0)
	  for i from 1 to (1- size)
	  for sum = (aref correlation 0 i)
	  and sumsq = (aref correlation i i)
	  for avg = (/ sum len)
	  for sigma = (sqrt (/  (- sumsq (* sum avg))
				len))
	  do
	     (assert (> sumsq (* avg sum)) () "Strange avg (~s = ~s / ~s) and sumsq (~s)" avg sum len sumsq)
	     (setf (aref A 0 i) (* -1s0 (/ avg sigma))
		   (aref A i i) (/ 1s0 sigma))
	  finally (return (values (times X a) A)))))
