(defparameter *matrix-field* 'single-float)
(defvar *default-declarations*
  '((scalar alpha beta))
  "Default list of declarations for with-matrixes.

If the head of an item is symbol SCALAR, interpret following symbols as scalars.
Otherwise it denotes name of a variable and follows its size.")

(defvar *matrix-zero* (coerce 0 *matrix-field*))

(defun compatible-size (s1 s2)
  (or (eq :any s1)
      (eq :any s2)
      (= s1 s2)))

(defmacro with-matrixes (expr
			 &key (declarations *default-declarations*)
			 (field *matrix-field*)
			 &environment env)
  "Run expr in linear algebra context.

Use optional declarations to indicate scalars or matrix sizes."
  (let ((*matrix-zero* (coerce 0 field))
	(*matrix-field* field)
	(*default-declarations* declarations))
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
		  (declare (optimize speed (debug 0) (safety 0))
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

(defun handle-summation (declarations expr env)
  (destructuring-bind (a-term . more-terms) expr
    (multiple-value-bind (a-res-type a-sizes a-res-expr a-assertions)
	(with-matrixes* declarations a-term env)
      (unless more-terms
	(return-from handle-summation
	  (values a-res-type a-sizes a-res-expr a-assertions)))
      (multiple-value-bind (b-res-type b-sizes b-res-expr b-assertions)
	  (with-matrixes* declarations `(+ ,@more-terms) env)
	(cond
	  ((and (eq 'scalar a-res-type)
		(eq 'scalar b-res-type))
	   (values 'scalar nil `(+ ,a-res-expr ,b-res-expr)
		   (append a-assertions b-assertions)))
	  ((or (eq 'scalar a-res-type) (eq 'scalar b-res-type))
	   (error "Cannot sum scalar and matrix: ~s ~s" a-term more-terms))
	  (t
	   (values 'matrix `(,(car a-sizes) ,(cadr a-sizes))
		   (lambda (i j)
		     `(+ ,(matrix-element-of a-res-type a-res-expr i j)
			 ,(matrix-element-of b-res-type b-res-expr i j)))
		   (append a-assertions b-assertions
			   `((assert (and
				      (compatible-size ,(car a-sizes) ,(car b-sizes))
				      (compatible-size ,(cadr a-sizes) ,(cadr b-sizes)))
				     ()
				     "Summation mismatch: ~s ~s" ',a-term ',more-terms))))))))))

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

(defun handle-diag (declarations expr env)
  (declare (ignore declarations env))
  (values 'matrix (list (length expr) (length expr))
	  (lambda (i j)
	    `(if (= ,i ,j) (elt ,(apply 'vector  expr) ,i) ,(coerce 0 *matrix-field*)))))

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

(defun handle-setq (declarations expr env)
  (destructuring-bind (place value) expr
      (multiple-value-bind (res-type res-size res-expr assertions)
	  (with-matrixes* declarations place env)
	(declare (ignore res-expr))
	(multiple-value-bind (src-type src-size src-expr src-assertions)
	    (with-matrixes* declarations value env)
	  (assert (symbolp place))
	  (cond
	    ((eq 'scalar res-type)
	     (assert (eq 'scalar src-type) () "Assigning non-scalar to scalar")
	     (values 'scalar nil `(setq ,(cadr expr) ,src-expr)
		     (append  src-assertions assertions)))
	    (t
	     (assert (not (eq 'scalar src-type)) ()
		     "Assigning a scalar ~s to a matrix ~s")
	     ;; returns 0
	     (values 'scalar res-size
		     (let ((i (gensym "I"))
			   (j (gensym "J")))
		       `(progn
			  (dotimes (,i ,(car res-size) res)
			    (declare (fixnum ,i))
			    (dotimes (,j ,(cadr res-size))
			      (declare (fixnum ,j))
			      (setf (aref ,place ,i ,j)
				    ,(matrix-element-of src-type src-expr i j))))
			  ,place)
		       (append assertions src-assertions
			       `((assert (and
					  (compatible-size ,(car src-size) ,(car res-size))
					  (compatible-size ,(cadr src-size) ,(cadr res-size)))
					 ()
					 "Cannot assign ~S to ~S"
					 ',value ',place)))))))))))

(defun with-matrixes* (declarations expr env)
  (when (atom expr)
    (return-from with-matrixes*
      (cond
	((arrayp expr) (values 'copy (array-dimensions expr) expr))
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
    (+ (handle-summation declarations (cdr expr) env))
    (trace (handle-trace declarations expr env))
    (transpose (handle-transpose declarations expr env))
    (setq (handle-setq declarations (cdr expr) env))
    (diag (handle-diag declarations (cdr expr) env))))

(defun matrix-element-of (type expr i j)
  (if (eq type 'copy) `(aref ,expr ,i ,j)
      (funcall expr i j)))
