(in-package sb-vm)

(defknown (f4+ f4- f4*) ((simd-pack single-float) (simd-pack single-float))
    (simd-pack single-float)
    (movable flushable always-translatable)
  :overwrite-fndb-silently t)

(define-vop (f4+)
  (:translate f4+)
  (:policy :fast-safe)
  (:args (x :scs (single-sse-reg) :target r) (y :scs (single-sse-reg)))
  (:arg-types simd-pack-single simd-pack-single)
  (:results (r :scs (single-sse-reg)))
  (:result-types simd-pack-single)
  (:generator 3
    (cond ((location= r y)
           (inst addps y x))
          (t
           (move r x)
           (inst addps r y)))))

(define-vop (f4*)
  (:translate f4*)
  (:policy :fast-safe)
  (:args (x :scs (single-sse-reg) :target r) (y :scs (single-sse-reg)))
  (:arg-types simd-pack-single simd-pack-single)
  (:results (r :scs (single-sse-reg)))
  (:result-types simd-pack-single)
  (:generator 4
    (cond ((location= r y)
           (inst mulps y x))
          (t
           (move r x)
           (inst mulps r y)))))

(macrolet ((stubize (a)
	     `(defun ,a (x y) (,a x y))))
  (stubize f4+)
  (stubize f4*))

(in-package linear-algebra)

(defun simd-scalar (a)
  (sb-vm::%make-simd-pack-single a a a a))

(defmacro with-simd-matrixes (op &rest args)
  `(regression::with-matrixes ,op ,@args
     :field (sb-vm::simd-pack single-float)
     :matrix-zero (simd-scalar 0.0)
     :adder sb-vm::f4+
     :multiplier sb-vm::f4*))

(locally (declaim (optimize (speed 3) (safety 0) (debug 0)))
  (defun test (a b)
    (declare ((simple-array (sb-vm::simd-pack single-float) (3 3)) a b))
    (with-simd-matrixes  (+ a a)
      :target a
      :use-assertions nil
      :optimize ((safety 0) (debug 0) (speed 3)))))

(let ((a (make-array '(3 3) :element-type '(sb-vm::simd-pack single-float)
			    :initial-element (sb-vm::%make-simd-pack-single 1.0 2.0 3.0 4.0))))
  (with-simd-matrixes  (+ a a)
    :optimize ((safety 0) (debug 0) (speed 3))))

(let (
      (a (make-array '(3 3) :element-type '(sb-vm::simd-pack single-float)
			   :initial-element (sb-vm::%make-simd-pack-single 1.0 2.0 3.0 4.0))))
  (with-simd-matrixes  (trace (* a (transpose a)))))
