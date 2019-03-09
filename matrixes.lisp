(defparameter *matrix-field* 'single-float)
(defvar *matrix-zero*)

(defmacro with-matrixes (declarations expr &environment env)
  (let ((*matrix-zero* (sb-pcl:class-prototype (find-class *matrix-field*))))
    (multiple-value-bind (res-type sizes res-expr assertions)
	(with-matrixes* declarations expr env)
      `(progn
	 ,@assertions
	 ,(ecase res-type
	    (copy res-expr)
	    (scalar res-expr)
	    (matrix
	     (let ((i (gensym "I"))
		   (j (gensym "J")))
	       `(let ((res (make-array (list ,@sizes) :element-type ',*matrix-field*
						      :initial-element ,*matrix-zero*)))
		  (declare (optimize speed)
			   ,@(mapcar (lambda (d)
				       `((simple-array ,*matrix-field* (,(cadr d) ,(caddr d)))
					 ,(car d)))
				     (remove 'scalar declarations :key 'car)))
		  (dotimes (,i ,(car sizes) res)
		    (declare (fixnum ,i))
		    (dotimes (,j ,(cadr sizes))
		      (declare (fixnum ,j))
		      (setf (aref res ,i ,j)
			    ,(funcall res-expr i j))))))))))))

(defun handle-multiply (declarations expr env)
  (multiple-value-bind (a-res-type a-sizes a-res-expr a-assertions)
      (with-matrixes* declarations (cadr expr) env)
    (unless (cddr expr)
      (return-from handle-multiply
	   (values a-res-type a-sizes a-res-expr a-assertions)))
       (multiple-value-bind (b-res-type b-sizes b-res-expr b-assertions)
	   (with-matrixes* declarations `(* ,@(cddr expr)) env)
	 (cond
	   ((and (eq 'scalar a-res-type)
		 (eq 'scalar b-res-type))
	    (values 'scalar nil `(* ,a-res-expr ,b-res-expr)
		    (append a-assertions b-assertions)))
	   ((eq 'scalar a-res-type)
	    (values 'matrix b-sizes
		    (lambda (i j)
		      (if (eq b-res-type 'copy) `(* ,a-res-expr (aref ,b-res-expr ,i ,j))
			  `(* ,a-res-expr ,(funcall  b-res-expr i j))))
		    (append a-assertions b-assertions)))
	   (t
	    (values 'matrix `(,(car a-sizes) ,(cadr b-sizes))
		    (lambda (i j)
		      (let ((k (gensym "K"))
			    (res (gensym "SUM")))
			`(let ((,res ,*matrix-zero*))
			   (declare (,*matrix-field* ,res))
			   (dotimes (,k ,(cadr a-sizes) ,res)
			     (incf ,res (* ,(matrix-element-of a-res-type a-res-expr i k)
					   ,(matrix-element-of b-res-type b-res-expr k j)))))))
		    (append a-assertions b-assertions
			    `((assert (= ,(cadr a-sizes) ,(car b-sizes))
				      ()
				      "Matrix multiply mismatch: ~sx~s" ',a-sizes ',b-sizes)))))))))


(defun handle-trace (declarations expr env)
  (multiple-value-bind (res-type sizes res-expr assertions)
      (with-matrixes* declarations (cadr expr) env)
    (assert (null (cddr expr)))
    (cond
      ((eq 'scalar res-type)
       (values 'scalar nil res-expr assertions))
      (t
       (values 'scalar nil
	       (let ((k (gensym "K"))
		     (res (gensym "SUM")))
		 `(let ((,res ,*matrix-zero*))
		    (declare (,*matrix-field* ,res))
		    (dotimes (,k ,(car sizes) ,res)
		      (incf ,res (* ,(matrix-element-of res-type res-expr k k))))))
	       (append assertions
		       `((assert (= ,(cadr sizes) ,(car sizes))
				 ()
				 "Trace on non-square: ~s" ',sizes ))))))))

(defun handle-transpose (declarations expr env)
  (multiple-value-bind (res-type sizes res-expr assertions)
      (with-matrixes* declarations (cadr expr) env)
    (assert (null (cddr expr)))
    (cond
      ((eq 'scalar res-type)
       (values 'scalar nil res-expr assertions))
      (t
       (values 'matrix (reverse sizes)
	       (lambda (i j)
		 (matrix-element-of res-type res-expr j i))
	       assertions)))))

(defun with-matrixes* (declarations expr env)
  (when (atom expr)
    (return-from with-matrixes*
      (cond
	((constantp expr env) (values 'scalar nil expr))
	((and (symbolp expr) (member expr (assoc 'scalar declarations)))
	 (values 'scalar nil expr))
	(t (values 'copy
		   (or (cdr (assoc expr declarations))
		       `((array-dimension ,expr 0)
			 (array-dimension ,expr 1)))
		   expr)))))
  (ecase (car expr)
    (* (handle-multiply declarations expr env))
    (trace (handle-trace declarations expr env))
    (transpose (handle-transpose declarations expr env))))

(defun matrix-element-of (type expr i j)
  (if (eq type 'copy) `(aref ,expr ,i ,j)
      (funcall expr i j)))

(with-matrixes nil 12)
(with-matrixes nil (* 12))
(with-matrixes nil (* 12 12))

(let ((B (make-array '(5 7) :element-type 'single-float :initial-element 3.0)))
  (with-matrixes ((B 5 7))
    (* 12.0 B)
    :field single-float))

(let ((B (make-array '(5 7) :element-type 'single-float :initial-element 3.0))
      (C (make-array '(7 5) :element-type 'single-float :initial-element 1.0))
      (p 0.1))
  (with-matrixes ((B 5 7) (C 7 5)
		  (scalar p))
    (trace (* p B C))))

(let ((B (make-array '(5 7) :element-type 'single-float :initial-element 3.0)))
  (with-matrixes ((B 5 7)
		  (scalar p))
    (* B (transpose B))))
