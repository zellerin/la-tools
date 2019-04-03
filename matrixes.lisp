(in-package linear-algebra)

(defparameter *matrix-field* 'single-float
  "Default field for the matrix elements.")

(defvar *matrix-optimize*
  (or nil '(speed (safety 0) (debug 0)))
  "Default optimization for the matrix calculation.")

(defvar *default-declarations*
  '((scalar alpha beta))
  "Default list of declarations for with-matrixes.

If the head of an item is symbol SCALAR, interpret following symbols as scalars.
Otherwise it denotes name of a variable and follows its size.")

(defvar *matrix-zero* (coerce 0 *matrix-field*)
  "Zero in the field of matrix elemets. It is rebound automatically inside `with-matrixes' macro.")

(define-condition matrix-error (simple-error)
  ((op     :accessor get-op     :initarg :op)
   (params :accessor get-params :initarg :params)
   (sizes  :accessor get-sizes  :initarg :sizes))
  (:documentation "Class of conditions that are raised for matrix calculations.")
  (:report
   (lambda (o stream)
     (format stream "Matrix error during ~A ~A" (get-op o) (get-params o)))))

(defun compatible-size (s1 s2)
  "Are the sizes of s1 and s2 compatible, that is, same or one or both is :any?"
  (or (eq :any s1) (eq :any s2) (= s1 s2)))

(defun assert-compatible-size (s1 s2 op params)
  "Test assertion that s1 and s2 are compatible before doing operation with parameters."
  (assert (compatible-size s1 s2)
	  ()
	  'matrix-error :op op :params params
	  :sizes (list s1 s2)))

(defmacro with-matrixes (expr
			 &key (declarations *default-declarations*)
			   (field *matrix-field*)
			   (optimize *matrix-optimize*)
			   target
			 &environment env)
  "Run expr in linear algebra context.

Use optional declarations to indicate scalars or matrix sizes."
  (let ((*matrix-zero* (coerce 0 field))
	(*matrix-field* field)
	(*default-declarations* declarations))
    (multiple-value-bind (res-type sizes res-expr assertions bindings)
	(with-matrixes* declarations expr env)
      `(let ,bindings
	 ,@assertions
	 ,(ecase res-type
	    (copy res-expr)
	    (scalar res-expr)
	    (matrix
	     (let ((i (gensym "I"))
		   (j (gensym "J")))
	       `(let ((res ,(or target
			       `(make-array (list ,@sizes) :element-type ',*matrix-field*
							  :initial-element ,*matrix-zero*))))
		  (declare (optimize ,@optimize)
			   ,@(mapcar (lambda (d)
				       `((simple-array ,*matrix-field* (,(cadr d) ,(caddr d)))
					 ,(car d)))
				     (remove-if (lambda (a) (member a '(vector scalar))) declarations
						:key 'car)))
		  (dotimes (,i ,(car sizes) res)
		    (dotimes (,j ,(cadr sizes))
		      (setf (aref res ,i ,j)
			    ,(funcall res-expr i j))))))))))))

(defun handle-summation (declarations expr env op)
  (cond
    ((eq (with-matrixes* declarations (car expr) env) 'scalar)
     (error "Summation of scalars not set yet"))
    (t (handle-map declarations `(,op ,@expr) env))))

(defun handle-multiply (declarations expr env)
  (multiple-value-bind (a-res-type a-sizes a-res-expr a-assertions a-bindings)
      (with-matrixes* declarations (car expr) env)
    (unless (cdr expr)
      (return-from handle-multiply
	   (values a-res-type a-sizes a-res-expr a-assertions a-bindings)))
       (multiple-value-bind (b-res-type b-sizes b-res-expr b-assertions b-bindings)
	   (with-matrixes* declarations `(* ,@(cdr expr)) env)
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
		    (append a-assertions b-assertions)
		    (append a-bindings b-bindings)))
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
				      'matrix-error :op "multiplication"
				      (list ',(car expr) ',(cdr expr)))))
		    (append a-bindings b-bindings)))))))

(defun handle-trace (declarations expr env)
  (multiple-value-bind (res-type sizes res-expr assertions bindings)
      (with-matrixes* declarations (car expr) env)
    (assert (null (cdr expr)))
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
				 "Trace on non-square: ~s" ',sizes )))
	       bindings)))))

(defun handle-map (declarations expr env)
  (multiple-value-bind (res-type sizes res-expr assertions bindings)
      (with-matrixes* declarations (cadr expr) env)
    (cond
      ((eq 'scalar res-type)
       (error 'matrix-error :op "map" :params expr))
      (t
       (loop for more-args in (cddr expr)
	     for (new-res-type new-sizes new-res-expr new-assertions new-bindings)
	       = (multiple-value-call 'list (with-matrixes* declarations more-args env))
	     append new-bindings into all-bindings
	     append (append new-assertions
			    `((assert-compatible-size ,(car new-sizes) ,(car sizes)
						      ',(car expr) ',(cdr expr))
			      (assert-compatible-size ,(cadr new-sizes) ,(cadr sizes)
						      ',(car expr) ',(cdr expr))))
	     into all-assertions
	     do
		(assert (not (equalp 'scalar new-res-type)))
	     collect new-res-expr into exprs
	     collect new-res-type into res-types
	     finally (return (values 'matrix
				     sizes
				     (lambda (i j)
				       `(,(car expr)
					 ,(matrix-element-of res-type res-expr i j)
					 ,@(mapcar (lambda (type expr)
						     (matrix-element-of type expr i j))
						   res-types exprs)))
				     (append assertions all-assertions)
				     (append bindings all-bindings))))))))

(defun handle-diag (declarations expr env)
  (declare (ignore declarations env))
  (values 'matrix (list (length expr) (length expr))
	  (lambda (i j)
	    `(if (= ,i ,j) (elt ,(apply 'vector  expr) ,i) ,(coerce 0 *matrix-field*)))))


(defun handle-transpose (declarations expr env)
  (multiple-value-bind (res-type sizes res-expr assertions bindings)
      (with-matrixes* declarations (car expr) env)
    (assert (null (cdr expr)))
    (cond
      ((eq 'scalar res-type)
       (values 'scalar nil res-expr assertions bindings))
      (t
       (values 'matrix (reverse sizes)
	       (lambda (i j)
		 (matrix-element-of res-type res-expr j i))
	       assertions bindings)))))

(defvar *with-matrixes-handlers*
  '((+ handle-summation +)
    (- handle-summation -)
    (map handle-map)
    (diag handle-diag)
    (trace handle-trace)
    (transpose handle-transpose)
    (* handle-multiply)))

(defmacro define-matrix-handler (symbol name &body code)
  `(progn
     (defun ,name (declarations expr env)
       ,@code)
     (pushnew (list ',symbol ',name) *with-matrixes-handlers* :key 'car)))

(define-matrix-handler aref handle-aref
  (values 'copy `((array-dimension (aref ,@expr) 0)
		  (array-dimension (aref ,@expr) 1))
	  `(aref ,@expr)))

(defun with-matrixes* (declarations expr env)
  (when (atom expr)
    (return-from with-matrixes*
      (cond
	((or
	  (vectorp expr)
	  (and (symbolp expr) (member expr (assoc 'vector declarations))))
	 (let ((rows (gensym (format nil "ROWS ~(~A~)" expr))))
	   (values 'matrix `(,rows 1)
		   (lambda (a b)
		     (declare (ignore b))
		     `(aref ,expr ,a))
		   nil
		   `((,rows (length ,expr))))))
	((arrayp expr) (values 'copy (array-dimensions expr) expr))
	((or (constantp expr env)
	     (and (symbolp expr) (member expr (assoc 'scalar declarations))))
	 (values 'scalar nil expr))
	(t (let ((rows (gensym (format nil "ROWS ~(~A~)" expr)))
		 (cols (gensym (format nil "COLS ~(~A~)" expr))))
	     (values 'copy
		     (or (cdr (assoc expr declarations))
			 `(,rows ,cols))
		     expr
		     nil
		     `((,rows (array-dimension ,expr 0))
		       (,cols  (array-dimension ,expr 1)))))))))
  (let ((handler (assoc (car expr) *with-matrixes-handlers*)))
    (assert handler () "No matrix handle for ~A" (car expr))
    (apply (cadr handler) declarations (cdr expr) env (cddr handler))))

(defun matrix-element-of (type expr i j)
  (if (eq type 'copy) `(aref ,expr ,i ,j)
      (funcall expr i j)))
