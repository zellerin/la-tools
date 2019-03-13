(in-package regression)

(declaim (inline sigma dsigma times-transposed))

(defun sigma (x)
  "Calculate logistic function x -> σ(x)

        1
σ = ---------
    1+exp(-x)

for parameter X. Care is taken not to underflow for
large negative values of x."
  (if (> x #.(- 1s0 (log most-positive-single-float)))
      (/ 1s0 (+ 1s0 (exp (- x))))
      0s0))

(defun dsigma (diff sigma)
  "Help calculate derivation of function of sigma, using formula

f(sigma(x))' = f'(sigma) . sigma . (1-sigma)."
  (* sigma (- 1s0 sigma) diff))

(defun trace-times-transposed (A B)
  "Return trace of matrix product, =Tr Aᵀ⋅B=

This is same as Sum(Aij Bij) over all i, j."
  (with-matrixes (trace (* (transpose A) B))))

(defun times-transposed (A B)
  "Return trace of matrix product, =Tr Aᵀ⋅B="
  (with-matrixes (* (transpose A) B)))

(defun make-random-array (x y scale)
  "Make random array of single floats in scale."
  (let ((res (make-array (list x y) :element-type 'single-float)))
    (dotimes (i x res)
      (dotimes (j y)
	(setf (aref res i j) (- scale (random (* 2 scale))))))))

(defun normalize (x)
  "Assuming that first column of X is all ones, return normalized matrix X and normalizer A, i.e.,
X -> (X' A) where
- X' = X.A,
- and average of each X column is 0 and standard deviation 1."
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
	  finally (return (values (with-matrixes (* X a)) A)))))
