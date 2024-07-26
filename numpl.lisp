(ql:quickload :petalisp)

(defpackage :numpl
  (:use :common-lisp :petalisp :petalisp.core)
  (:local-nicknames (:pl :petalisp)))

(in-package :numpl)


;; * Array creation
(defun zeros (shape)
  "Return lazy-array full of zeros with shape SHAPE"
  (pl:lazy-reshape 0.0 shape))

(defun ones (shape)
  "Return lazy-array full of ones with shape SHAPE"
  (pl:lazy-reshape 1.0 shape))

(defun full (shape fill-value)
  "Return lazy-array of shape SHAPE filled with FILL-VALUE"
  (pl:lazy-reshape fill-value shape))

(defun cyclically-shift (arr num axis)
  "Cyclically shift elements of array ARR along AXIS NUM times"
  (pl:with-lazy-arrays (arr)
    (when (> axis (1- (length (pl:lazy-array-dimensions arr))))
      (error "Invalid axis ~S for shape ~S." axis (pl:lazy-array-dimensions arr)))
    (let* ((dims (pl:lazy-array-dimensions arr))
           (axis-len (nth axis dims))
           (split-point (mod  (- axis-len num) axis-len))
           (end-component (pl:lazy-slices arr (range split-point axis-len) axis)))
      (pl:lazy-stack axis end-component (pl:lazy-slices arr (range 0 split-point) axis)))))


;;;TODO
(defun swap-axes (arr axis1 axis2))


;;; TODO write this function
(defun eye (N &optional (M nil) (k 0))
  "Return 2D lazy-array with ones on diagnoal, zeros elsewhere"
  (when (null M)
    (setf M N))
  (let ((max-dim (max N M)))))


(defun zeros-like (arr)
  "Return lazy-array full of zeros with same shape as ARR"
  (cond ((arrayp arr) (pl:lazy-reshape 0.0 (pl:array-shape arr)))
        ((pl:lazy-array-p arr) (pl:lazy-reshape 0.0 (pl:lazy-array-shape arr)))))

(defun ones-like (arr)
  "Return lazy-array full of ones with same shape as ARR"
  (cond ((arrayp arr) (pl:lazy-reshape 1.0 (pl:array-shape arr)))
        ((pl:lazy-array-p arr) (pl:lazy-reshape 1.0 (pl:lazy-array-shape arr)))))


(defun full-like (arr fill-value)
  "Return lazy-array full of FILL-VALUE with same shape as ARR"
  (cond ((arrayp arr) (pl:lazy-reshape fill-value (pl:array-shape arr)))
        ((pl:lazy-array-p arr) (pl:lazy-reshape fill-value (pl:lazy-array-shape arr)))))
