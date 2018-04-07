
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

(defparameter *multipliers* nil
  "List of discriminating functions and associated multipliers.")

(defparameter *linear-combinations* nil
  "List of discriminating functions and associated x=ax+by functions.")

(defparameter *transpose-multipliers* nil
  "List of discriminating functions and associated x=ax+by functions.")

(macrolet
    ((def (element zero a-type &key (b-type a-type) (c-type a-type) (speed 3))
       `(progn
	  ,@(loop for a-order in '((row item) (item row))
		  and b-order in '((item col) (item col))
		  and target in '(*multipliers* *transpose-multipliers*)
		  collect
		  `(push (cons
			  (lambda (a b res)
			    (and (typep a ',a-type)
				 (typep b ',b-type)
				 (typep res ',c-type)))
			  (lambda (A B res rows batch-size cols)
			    (declare (optimize (speed ,speed) (safety 2) (debug 0))
				     (,a-type A)
				     (,b-type B)
				     (,c-type res)
				     (fixnum rows batch-size cols))
			    (dotimes (col cols res)
			      (dotimes (row rows)
				(let ((val ,zero))
				  (declare (,element val))
				  (dotimes (item batch-size)
				    (incf val
					  (* (aref A ,@a-order)
					     (aref B ,@b-order)))
				    (setf (aref res row col) val)))))))
			 ,target)))))

  (def t 0 (array) :speed 1)			; base case
  (def single-float 0s0 (simple-array single-float))
  (def single-float 0s0 (simple-array single-float (500 500)))
  (def double-float 0d0 (simple-array double-float)))

(macrolet
    ((def (element v-type &key (speed 3))
       `(push (cons
	       (lambda (res x y a b)
		 (and (typep x ',v-type)
		      (typep y ',v-type)
		      (typep res ',v-type)
		      (typep a ',element)
		      (typep b ',element)))
	       (lambda (res a x b y rows cols)
		 (declare (optimize (speed ,speed) (safety 1) (debug 0))
			  (,v-type x y res)
			  (,element a b)
			  (fixnum rows cols))
		 (dotimes (col cols res)
		   (dotimes (row rows)
		     (setf (aref res row col)
			   (+ (* a (aref x row col))
			      (* b (aref y row col))))))))
	      *linear-combinations*)))

  (def t (array) :speed 1) ; base case
  (def single-float (simple-array single-float))
  (def single-float (simple-array single-float (500 500)))
  (def double-float (simple-array double-float)))

(defun times-into (A B res cols batch-size rows)
  (loop for (check . fn) in *multipliers*
	when (funcall check a b res)
	  do (return (funcall fn a b res cols batch-size rows))
	finally (error "No matching fn")))

(defun times-transposed-into (A B res cols batch-size rows)
  (loop for (check . fn) in *transpose-multipliers*
	when (funcall check a b res)
	  do (return (funcall fn a b res cols batch-size rows))
	finally (error "No matching fn")))

(defun linear-combination-into (res a x b y)
  (assert (equalp (array-dimensions x) (array-dimensions y)))
  (loop for (check . fn) in *linear-combinations*
	when (funcall check res x y a b)
	  do (return (funcall fn res a x b y
			      (array-dimension x 0)
			      (array-dimension x 1)))
	finally (error "No matching fn")))

(defun times (A B)
  (let* ((batch-size (array-dimension A 1))
	 (rows (array-dimension A 0))
	 (cols (array-dimension B 1)))
    (assert (= batch-size (array-dimension B 0)))
    (times-into A B (make-array (list rows cols)
				:element-type (array-element-type A))
		rows batch-size cols)))

(defun times-transposed (A B)
  (let* ((batch-size (array-dimension A 0))
	 (cols (array-dimension A 1))
	 (rows (array-dimension B 1)))
    (assert (= batch-size (array-dimension B 0)))
    (times-transposed-into A B (make-array (list cols rows)
				:element-type (array-element-type A))
			   cols batch-size rows)))

(defun linear-update (a x b y)
  (linear-combination-into x a x b y))

(defun linear-combination (a x b y)
  (linear-combination-into
   (make-array (array-dimensions x)
	       :element-type (array-element-type x))
   a x b y))

