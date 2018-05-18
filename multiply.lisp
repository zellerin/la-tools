
(in-package linear-algebra)

;;;; Multiplication functions on matrixes of three variants:
;;;; - AB (i.e., A(i j) B(j i)
;;;; - ABt
;;;;
;;;; Each list contains several possible funtions to multiply, each
;;;; specialing on some matrix types.
;;;;
;;;; Possible extensions: allow multiply vector and matrix in order,
;;;; apply function in order, etc. - but you usually know that you
;;;; want to do this, so you can specialize yourself.
;;;;
;;;; Typically, it make sense to majke a generic function (without
;;;; speed optimization) and on simple-array of single-float or double-float
;;;;
;;;; You can get even faster if you specify array size; this is not
;;;; tested only for square matrixes.
;;;;
;;;; Some speed optimizations are in place, some not. In particular,
;;;; no effort to optimize memory access was done, as experiments did
;;;; not show it helpful on my machine.


(defun update-matrix (matrix fn rows cols)
  "Create a matrix array based on fn(row, col)"
  (dotimes (col cols matrix)
    (dotimes (row rows)
      (setf (aref matrix row col)
	    (funcall fn row col)))))

(defun trace-by-fn (fn rows)
  (let ((res 0))
    (dotimes (i rows res)
      (incf res (funcall fn i i)))))

(defun make-fn-list (a-order b-order)
  (flet
      ((compute-specialized (element zero
		 &key (a-type `(simple-array ,element (* *)))
		   (b-type a-type) (speed 3))
	 (cons
	  (lambda (a b &rest more)
	    (declare (ignore more))
	    (and (typep a a-type)
		 (typep b b-type)))
	  (compile nil
		   `(lambda (A B batch-size)
		      (declare (optimize (speed ,speed)
					 (safety 2)
					 (debug 0))
			       (,a-type A)
			       (,b-type B)
			       (fixnum batch-size))

		      (flet ((compute-cell (row col)
			       (let ((val ,zero))
				 (declare (,element val))
				 (dotimes (item batch-size val)
				   (incf val
					 (* (aref A ,@a-order)
					    (aref B ,@b-order)))))))
			;; inline helps with known matix size, but in general not.
			#'compute-cell
			#+nil(update-matrix res #'compute-cell rows cols)))))))
    (mapcar (lambda (a) (apply #'compute-specialized a))
	    '((single-float 0s0 :a-type (simple-array single-float (500 500)))
	      (single-float 0s0)
	      (t 0 :a-type array :speed 1)))))

(defun make-linear-combination ()
  (flet ((def (element &key (v-type `(simple-array ,element (* *))) (speed 3))
	   (cons
	    (lambda (res a x b y &rest more)
	      (declare (ignore more))
	      (and (typep x v-type)
		   (typep y v-type)
		   (typep res v-type)
		   (typep a element)
		   (typep b element)))
	    (compile nil
		     `(lambda (res a x b y rows cols)
		       (declare (optimize (speed ,speed) (safety 1) (debug 0))
				(,v-type x y res)
				(,element a b)
				(fixnum rows cols))
			   (update-matrix res
					  (lambda (row col)
					    (+ (* a (aref x row col))
					       (* b (aref y row col))))
					  rows cols))))))

    (list
     (def 'single-float :v-type '(simple-array single-float (500 500)))
     (def 'single-float)
     (def t :v-type '(array) :speed 1))))

(defun find-applicable-fn (pars candidates)
   (cdr (or
	 (find-if (lambda (a) (apply a pars))
		  candidates :key 'car)
	 (error "No matching function")))  )

(defun call-applicable-fn (pars candidates)
  (apply
   (find-applicable-fn pars candidates)
   pars))

(defun times-into (rows cols res &rest pars)
  (update-matrix res
		 (call-applicable-fn pars
				     (load-time-value
				      (make-fn-list '(row item) '(item col))))
		 rows cols))

(let ((fns (load-time-value (make-fn-list '(item row) '(item col)))))
  (defun times-transposed-into (rows cols res &rest pars)
    (update-matrix res
		   (call-applicable-fn pars fns)
		   rows cols))

  (defun trace-times-transposed* (rows &rest pars)
    (trace-by-fn (call-applicable-fn pars fns)
		 rows)))

(defun times-rev-transposed-into (rows cols res &rest pars)
  (update-matrix res
		 (call-applicable-fn pars
				     (load-time-value (make-fn-list '(row item) '(col item))))
		 rows cols))

(defun linear-combination-into (&rest pars)
  (call-applicable-fn pars
		      (load-time-value  (make-linear-combination))))

(defun times (A B)
  "Return matrix product =A⋅B=. Nondestructive.

Matrixes are represented as arrays both of input and output.
Signal error if A and B are not matrices, or if the number of rows in A does not match columns in B."
  (let* ((batch-size (array-dimension A 1))
	 (rows (array-dimension A 0))
	 (cols (array-dimension B 1)))
    (assert (= batch-size (array-dimension B 0)))
    (times-into rows cols
		(make-array (list rows cols)
			    :element-type (array-element-type A))
		A B 
		batch-size)))

(defun times-rev-transposed (A B)
  "Return matrix product =A⋅Bᵀ=. Nondestructive."

  (let* ((batch-size (array-dimension A 1))
	 (rows (array-dimension A 0))
	 (cols (array-dimension B 0)))
    (assert (= batch-size (array-dimension B 1)))
    (times-rev-transposed-into rows cols
			        (make-array (list rows cols)
					    :element-type (array-element-type A))
				A B
				batch-size)))

(defun times-transposed (A B)
  "Return matrix product =Aᵀ⋅B=. Nondestructive."
  (let* ((batch-size (array-dimension A 0))
	 (rows (array-dimension A 1))
	 (cols (array-dimension B 1)))
    (assert (= batch-size (array-dimension B 0)))
    (times-transposed-into rows cols
			   (make-array (list rows cols)
				       :element-type (array-element-type A))
			   A B
			   batch-size)))

(defun trace-times-transposed (A B)
  "Return trace of matrix product, =Tr Aᵀ⋅B="
  (let* ((batch-size (array-dimension A 0))
	 (rows (array-dimension A 1))
	 (cols (array-dimension B 1)))
    (assert (= batch-size (array-dimension B 0)))
    (assert (= cols rows))
    (trace-times-transposed* rows
			     A B
			     batch-size)))

(defun linear-update (a x b y)
  "Modify X to hold =aX+bY= and return it."
  (assert (equalp (array-dimensions x) (array-dimensions y)))
  (linear-combination-into x a x b y (array-dimension x 0)
			   (array-dimension x 1)))

(defun linear-combination (a x b y)
  "Return =aX+bY=. Nondestructive."
  (assert (equalp (array-dimensions x) (array-dimensions y)))
  (linear-combination-into
   (make-array (array-dimensions x)
	       :element-type (array-element-type x))
   a x b y (array-dimension x 0) (array-dimension x 1)))

(defun M- (a b)
  "Return matrix difference, =A-B="
  (linear-combination 1s0 a -1s0 b))
