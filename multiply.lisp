(defparameter *multipliers* nil
  "List of discriminating functions and associated multipliers.")

(defparameter *linear-combinations* nil
  "List of discriminating functions and associated x=ax+by functions.")

(defparameter *transpose-multipliers* nil
  "List of discriminating functions and associated x=ax+by functions.")

(defparameter *transpose-b-multipliers* nil
  "List of discriminating functions and associated x=ax+by functions.")

(macrolet
    ((def (element zero a-type &key (b-type a-type) (c-type a-type) (speed 3))
       `(progn
	  ,@(loop for a-order in '((col item) (item col) (col item))
		  and b-order in '((item row) (item row) (row item))
		  and target in '(*multipliers* *transpose-multipliers*
				  *transpose-b-multipliers*)
		  collect
		  `(push (cons
			  (lambda (a b res)
			    (and (typep a ',a-type)
				 (typep b ',b-type)
				 (typep res ',c-type)))
			  (lambda (A B res cols batch-size rows)
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
				    (setf (aref res col row) val)))))))
			 ,target)))))

  (def t 0 (array) :speed 1)			; base case
  (def single-float 0s0 (simple-array single-float))
  (def single-float 0s0 (simple-array single-float (500 500)))
  (def double-float 0d0 (simple-array double-float)))

(macrolet
    ((def (element v-type &key (speed 3))
       `(push (cons
	       (lambda (x y a b)
		 (and (typep x ',v-type)
		      (typep y ',v-type)
		      (typep a ',element)
		      (typep b ',element)))
	       (lambda (a x b y rows cols)
		 (declare (optimize (speed ,speed) (safety 1) (debug 0))
			  (,v-type x y)
			  (,element a b)
			  (fixnum rows cols))
		 (dotimes (col cols x)
		   (dotimes (row rows)
		     (setf (aref x row col)
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

(defun times-transposed-b-into (A B res cols batch-size rows)
  (loop for (check . fn) in *transpose-b-multipliers*
	when (funcall check a b res)
	  do (return (funcall fn a b res cols batch-size rows))
	finally (error "No matching fn")))

(defun linear-combination (a x b y)
  (assert (= (array-dimension x 0) (array-dimension y 0)))
  (assert (= (array-dimension x 1) (array-dimension y 1)))
  (loop for (check . fn) in *linear-combinations*
	when (funcall check x y a b)
	  do (return (funcall fn a x b y
			      (array-dimension x 0)
			      (array-dimension x 1)))
	finally (error "No matching fn")))

(defun times (A B)
  (let* ((batch-size (array-dimension A 1))
	 (cols (array-dimension A 0))
	 (rows (array-dimension B 1)))
    (assert (= batch-size (array-dimension B 0)))
    (times-into A B (make-array (list cols rows)
				:element-type (array-element-type A))
		cols batch-size rows)))

(defun times-transposed (A B)
  (let* ((batch-size (array-dimension A 0))
	 (cols (array-dimension A 1))
	 (rows (array-dimension B 1)))
    (assert (= batch-size (array-dimension B 0)))
    (times-transposed-into A B (make-array (list cols rows)
				:element-type (array-element-type A))
			   cols batch-size rows)))

(defun times-transposed-b (A B)
  (let* ((batch-size (array-dimension A 1))
	 (cols (array-dimension A 0))
	 (rows (array-dimension B 0)))
    (assert (= batch-size (array-dimension B 1)))
    (times-transposed-b-into A B (make-array (list cols rows)
				:element-type (array-element-type A))
		cols batch-size rows)))
