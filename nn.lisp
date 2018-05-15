(in-package regression)

;; Neural network is represented by a vector of transformation matrixes.
;; Last transformation is logistic or linear, other are logistic

(define-pair nn-estimate (X0 As)
    "Propagate X values across neural network with last layer ~(~A~).

Return array of X_i+1 and final Y"
  (loop
    with len = (length As)
    with res = (make-array (1+ len))
    initially (setf (aref res 0) X0)
    for i from 0 to len
    for A across As
    for X = (aref res i)
    do (setf (aref res (1+ i))
	     (if (= i len)
		 (estimate X A)
		 (logistic-estimate X A)))
    finally (return res)))

(defun linear-propagate-error (err A Y)
  (declare (ignore Y))
  (times-rev-transposed A err))

(defun logistic-propagate-error (err A Y)
  (times-rev-transposed
   (apply-fn2 #'float-dsigma (copy-array err) Y) A))

(eval-when (:compile-toplevel)
  (pushnew  'propagate-error *pairs*))

(define-pair nn-backpropagate (Yn Xs As sigmas rho)
    "Backpropagate neural network with last term ~(~A~)"
    (loop with len = (length As)
	  with orig-err  = (M- (aref Xs len) Yn)
	  with err = orig-err
	  for i from len downto 1
	  for sigma across sigmas
	  for A = (aref As (1- i))
	  for Y = (aref Xs i)
	  for X = (aref Xs (1- i))
	  for grad-A = (grad-A X (copy-array err) Y)
	  do
	     (linear-update rho A sigma grad-A)
	     (setq err
		   (if (= i (1- len))
			       (propagate-error err A Y)
			       (logistic-propagate-error err A Y)))

	  finally (return (times-transposed orig-err orig-err))))

;;; Test case
(defun test-nn ()
  (let* ((X (make-random-array 10 5 1s0))
	 (Y (aref
	     (linear-nn-estimate X
				 (vector (make-random-array 5 2 1s0)
					 (make-random-array 2 1 1s0)))
	     2))
	 (A (vector (make-random-array 5 2 1s0)
		    (make-random-array 2 1 1s0))))
    (loop for i from 0 to 3000
	  for Xs = (linear-nn-estimate X A)
	  do (print (linear-nn-backpropagate Y Xs A #(-1s0 -1s0) 0.9999)))))
