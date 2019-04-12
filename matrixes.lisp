(in-package linear-algebra)

(defparameter *matrix-field* 'single-float
  "Default field for the matrix elements.")

(defvar *matrix-optimize*
  (and nil '(speed (safety 1) (debug 0)))
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

(defclass base-matrix-object ()
  ()
  (:documentation "Mixin for expression objects that represent a matrix (as oposed to, e.g. scalar)"))

(defclass expr-with-children ()
  ((children   :accessor get-children   :initarg :children)
   (assertions :accessor get-assertions :initarg :assertions))
  (:default-initargs :children nil :assertions nil)
  (:documentation "Mixin for expression objects that depends on other objects (aka children).

Binding for chidren are prepended before bindings for the object, and
there might be assertions to ensure that the children are
compatible."))

(defclass expr-object (expr-with-children base-matrix-object)
  ((sizes      :accessor get-sizes      :initarg :sizes)
   (expr       :accessor get-expr       :initarg :expr))
  (:default-initargs)
  (:documentation "Object that represents a matrix with element values that need to be calculated.

The calculation is delayed until assigning to a target matrix or similar time."))

(defclass scalar-expr-object (expr-with-children)
  ((type       :accessor get-type       :initarg :type :allocation :class)
   (expr       :accessor get-expr       :initarg :expr)
   (sizes      :accessor get-sizes      :initarg :sizes :allocation :class))
  (:default-initargs :sizes nil)
  (:documentation "Object that represents a scalar expression (scalar value, trace of matrix, etc)"))

(defclass matrix-copyable-object (base-matrix-object)
  ((object       :accessor get-object       :initarg :object :reader get-expr)
   (row-var-name :accessor get-row-var-name :initarg :row-var-name)
   (col-var-name :accessor get-col-var-name :initarg :col-var-name))
  (:default-initargs :row-var-name (gensym "ROW") :col-var-name (gensym "COL"))
  (:documentation "Object that represents an expression that produces an array.

The size values of elements is determined at run time."))

(defclass matrix-literal-object (base-matrix-object)
  ((object       :accessor get-object       :initarg :object :reader get-expr))
  (:documentation "Object that represents a literal array (with known size and elements"))

(defgeneric matrix-element-of* (object i j)
  (:method ((object matrix-copyable-object) i j)
    `(aref ,(get-object object) ,i ,j))
  (:method ((object matrix-literal-object) i j)
    `(aref ,(get-object object) ,i ,j))
  (:method ((object expr-object) i j)
    (funcall (get-expr object) i j)))


(defun make-expr-object (sizes expr &optional assertions children)
  (make-instance 'expr-object
		     :sizes sizes
		     :expr expr
		     :assertions assertions
		     :children children))

(defgeneric get-bindings (object)
  (:method (object) nil)
  (:method ((object matrix-copyable-object))
    (let ((expr (get-object object)))
      `((,(get-row-var-name object) (array-dimension ,expr 0))
	(,(get-col-var-name object) (array-dimension ,expr 1)))))
  (:method :around ((object expr-with-children))
    (append (reduce 'append (get-children object) :key 'get-bindings)
	    (call-next-method))))

(defmethod get-sizes ((object matrix-copyable-object))
  (list (get-row-var-name object) (get-col-var-name object)))

(defmethod get-sizes ((object matrix-literal-object))
  (array-dimensions (get-object object)))

(defmethod get-assertions ((object matrix-literal-object))
  nil)

(defgeneric get-assertions (object)
  (:method (object) nil)
  (:method :around ((object expr-with-children))
    (append (reduce 'append (get-children object) :key 'get-assertions)
	    (call-next-method))))

(defgeneric assign-expr (object target)
  (:method (object (target symbol)) `(setf ,target ,(get-expr object)))
  (:method (object (target (eql nil))) (get-expr object))
  (:method ((object expr-object) res)
    (let ((i (gensym "I"))
	  (j (gensym "J")))
      `(dotimes (,i ,(car (get-sizes object)) ,res)
	(dotimes (,j ,(cadr (get-sizes object)))
	  (setf (aref ,res ,i ,j)
		,(funcall (get-expr object) i j))))))
  (:method ((object expr-object) (target (eql nil)))
    `(let ((target (make-array (list ,@(get-sizes object))
				:element-type ',*matrix-field*
				:initial-element ,*matrix-zero*)))
       ,(assign-expr object 'target))))

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
    (let ((result (calculate-matrix-object declarations expr env)))
      `(let ,(get-bindings result)
	 (declare (optimize ,@optimize)
		  ,@(mapcar (lambda (d)
			      `((simple-array ,*matrix-field* (,(cadr d) ,(caddr d)))
				,(car d)))
			    (remove-if (lambda (a) (member a '(vector scalar))) declarations
				       :key 'car)))
	 ,@(get-assertions result)

	 ,(assign-expr result target)))))

(defgeneric handle-multiply* (a b)
  (:method (a b) nil)
  (:method (a (b (eql nil))) a)
  (:method ((a scalar-expr-object) (b scalar-expr-object))
    (make-instance 'scalar-expr-object
		   :expr `(* ,(get-expr a) ,(get-expr b))
		   :children (list a b)))
  (:method ((a scalar-expr-object) (b base-matrix-object))
    (make-expr-object (get-sizes b)
		      (lambda (i j)
			`(* ,(get-expr a)
			    ,(matrix-element-of* b i j)))
		      nil (list a b)))
  (:method ((a base-matrix-object) (b base-matrix-object))
    (let ((a-sizes (get-sizes a))
	  (b-sizes (get-sizes b)))
      (make-expr-object `(,(car a-sizes) ,(cadr b-sizes))
			;; todo:define do-matrix-elements to create
			;; the per-matrix loop; then, sparse matrixes
			;; are easier
			(lambda (i j)
			  (let ((k (gensym "K"))
				(res (gensym "SUM")))
			    `(let ((,res ,*matrix-zero*))
			       (declare (,*matrix-field* ,res))
			       (dotimes (,k ,(cadr a-sizes) ,res)
				 (incf ,res (* ,(matrix-element-of* a i k)
					       ,(matrix-element-of* b k j)))))))
			`((assert (= ,(cadr a-sizes) ,(car b-sizes))
				  ()
				  'matrix-error :op "multiplication"
				  ))
			(list a b))))
  ;; optimization that can be removed: compute during compilation phase
  (:method ((a matrix-literal-object) (b matrix-literal-object))
    (make-instance 'matrix-literal-object
		   :object (funcall (compile ()
					     `(lambda () ,(assign-expr (call-next-method) nil)))))))

(defun handle-multiply (declarations expr env)
  (let ((a-obj (calculate-matrix-object declarations (car expr) env))
	(b-obj (and (cdr expr) (calculate-matrix-object declarations `(* ,@(cdr expr)) env))))
    (handle-multiply* a-obj b-obj)))

(defun handle-trace (declarations expr env)
  (assert (null (cdr expr)))
  (let* ((par-object (calculate-matrix-object declarations (car expr) env))
	 (sizes (get-sizes par-object)))
    (etypecase par-object
      (scalar-expr-object par-object)
      (base-matrix-object
       (make-instance 'scalar-expr-object
		    :expr (let ((k (gensym "K"))
				(res (gensym "SUM")))
			    `(let ((,res ,*matrix-zero*))
			       (declare (,*matrix-field* ,res))
			       (dotimes (,k ,(car sizes) ,res)
				 (incf ,res (* ,(matrix-element-of* par-object k k))))))
		    :assertions `((assert (= ,(cadr sizes) ,(car sizes))
					  ()
					  "Trace on non-square: ~s" ',sizes ))
		    :children (list par-object))))))

(defvar *testob*)

(defun handle-map (declarations expr env)
  (let* ((first-par
	   (calculate-matrix-object declarations (cadr expr) env))
	 (sizes (get-sizes first-par)))
    (etypecase first-par
      (scalar-expr-object
       (error 'matrix-error :op "map" :params expr))
      (base-matrix-object
       (loop with more-args-objs = (mapcar (lambda (a)
					     (calculate-matrix-object declarations a env))
					   (cddr expr))
	     for more-args in more-args-objs
	     append
	     `((assert-compatible-size ,(car (get-sizes  more-args))
				       ,(car sizes)
				       ',(car expr) ',(cdr expr))
	       (assert-compatible-size ,(cadr (get-sizes more-args)) ,(cadr sizes)
				       ',(car expr) ',(cdr expr)))
	       into all-assertions
	     do
		(assert (typep more-args 'base-matrix-object))
	     finally (return (make-expr-object sizes
					       (lambda (i j)
						 `(,(car expr)
						   ,(matrix-element-of* first-par i j)
						   ,@(mapcar (lambda (obj)
							       (matrix-element-of* obj i j))
							     more-args-objs)))
					       all-assertions
					       (cons first-par more-args-objs))))))))

(defvar *with-matrixes-handlers*
  '((map handle-map)
    (trace handle-trace)
    (transpose handle-transpose)
    (* handle-multiply)))

(defmacro define-matrix-handler (symbol name &body code)
  `(progn
     (defun ,name (declarations expr env)
       ,@code)
     (pushnew (list ',symbol ',name) *with-matrixes-handlers* :key 'car)))

(define-matrix-handler aref handle-aref
  (declare (ignore env declarations))
  (make-instance 'matrix-copyable-object
		 :object `(aref ,@expr)))

(define-matrix-handler + handle-plus
  (handle-map declarations `(+ ,@expr) env))

(define-matrix-handler - handle-minus
  (handle-map declarations `(- ,@expr) env))

(define-matrix-handler transpose handle-transpose
  (assert (null (cdr expr)))
  (let ((result (calculate-matrix-object declarations (car expr) env)))
    (etypecase result
      (scalar-expr-object expr)
      (base-matrix-object
       (make-expr-object (reverse (get-sizes result))
			 (lambda (i j)
			   (matrix-element-of* result j i))
			 nil (list result))))))



(defun calculate-matrix-object (declarations expr env)
  (when (atom expr)
    (return-from calculate-matrix-object
      (cond
	((or
	  (vectorp expr)
	  (and (symbolp expr) (member expr (assoc 'vector declarations))))
	 (let ((rows (gensym (format nil "ROWS ~(~A~)" expr))))
	   (make-expr-object `(,rows 1)
		   (lambda (a b)
		     (declare (ignore b))
		     `(aref ,expr ,a))
		   nil
		   `((,rows (length ,expr))))))
	((arrayp expr)
	 (make-instance 'matrix-literal-object :object expr))
	((or (constantp expr env)
	     (and (symbolp expr) (member expr (assoc 'scalar declarations))))
	 (make-instance 'scalar-expr-object :expr expr))
	(t
	 (make-instance 'matrix-copyable-object :object expr)))))
  (let ((handler (assoc (car expr) *with-matrixes-handlers*)))
    (assert handler () "No matrix handle for ~A" (car expr))
    (apply (cadr handler) declarations (cdr expr) env (cddr handler))))
