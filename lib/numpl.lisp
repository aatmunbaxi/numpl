(defpackage :numpl
  (:use
   :common-lisp
   :alexandria
   :petalisp)
  (:local-nicknames (nil :petalisp)
                    (:ax :alexandria))
  (:export
   ;; Array properties
   #:ndim
   #:shape
   #:size
   ;; Array creation
   #:zeros
   #:ones
   #:full
   #:random-array
   #:zeros-like
   #:ones-like
   #:full-like
   #:array-coordinate-grid
   #:identity-matrix
   #:arange
   #:linspace
   #:logspace
   ;; Array manipulation
   #:lazy-concat
   #:shift-entries
   #:cyclically-shift
   #:permute-dims
   #:transpose
   #:matrix-transpose
   #:swap-axes
   #:flip))

(in-package :numpl)



;;; Array properties
(defun ndim (arr)
  "Return number of dimensions for ARR"
  (with-lazy-arrays (arr)
    (lazy-array-rank arr)))

(defun shape (arr)
  "Return shape of ARR."
  (with-lazy-arrays (arr)
    (lazy-array-dimensions arr)))

(defun size (arr &optional (axis 0))
  "Return number of elements along a given axis"
  (with-lazy-arrays (arr)
    (nth axis (lazy-array-dimensions arr))))




;;; Array creation routines
;;;; Prefilled arrays
(defun zeros (shape &optional (element-type 'integer))
  "Return lazy-array full of zeros with shape SHAPE"
  (lazy-reshape (coerce 0 element-type) shape))

(defun ones (shape &optional (element-type 'integer))
  "Return lazy-array full of ones with shape SHAPE"
  (lazy-reshape (coerce 1 element-type) shape))

(defun full (shape fill-value &optional (element-type nil))
  "Return lazy-array of shape SHAPE filled with FILL-VALUE"
  (unless element-type
    (setf element-type (type-of fill-value)))
  (lazy-reshape (coerce fill-value element-type) shape))



(defun random-array (shape &optional (limit 256) (element-type nil))
  "Return lazy-array with SHAPE full of random values
with limit LIMIT."
  (let* ((element-type (unless element-type (typecase (type-of limit)
                                              (list (car (type-of limit)))
                                              (symbol (type-of limit)))))
         (array (make-array (shape-dimensions shape) )))
    (loop for index below (array-total-size array) do
      (setf (row-major-aref array index)
            (random limit)))
    (lazy (lambda (x) (coerce x element-type)) array)))


;; TODO write eye function


(defun zeros-like (arr &optional (element-type nil))
  "Return lazy-array full of zeros with same shape as ARR."
  (with-lazy-arrays (arr)
    (unless element-type
      (setf element-type (lazy-array-element-type arr)))
    (lazy-reshape (coerce 0 element-type) (lazy-array-shape arr))))


(defun ones-like (arr &optional (element-type nil))
  "Return lazy-array full of ones with same shape as ARR."
  (with-lazy-arrays (arr)
    (unless element-type
      (setf element-type (lazy-array-element-type arr)))
    (lazy-reshape (coerce 1 element-type) (lazy-array-shape arr))))

(defun full-like (arr fill-value &optional (element-type nil))
  "Return lazy-array full of FILL-VALUE with same shape as ARR."
  (with-lazy-arrays (arr)
    (unless element-type
      (setf element-type (lazy-array-element-type arr)))
    (lazy-reshape (coerce fill-value element-type) (lazy-array-shape arr))))


(defun array-coordinate-grid (shape)
  "Return lazy array whose elements are lists (i j k .. ) where i, j, k, ... etc.
are the coordinates at that point of the array with shape SHAPE."
  (apply #'lazy #'list (loop for axis
                               below (shape-rank shape)
                             collect
                             (lazy-index-components shape axis))))


(defun identity-matrix (N)
  (apply #'lazy-stack 0 (loop for i below N
                              collect
                              (lazy-overwrite (zeros (~ i (1+ i) ~ N))
                                              (lazy-reshape 1 (~ i (1+ i) ~ i (1+ i)))))))

;;;;; Arrays with ranges
(defun arange (stop &key (start 0)  (num 50))
  "Return lazy endpoint-exclusive, evenly-spaced range of numbers, starting from START
and ending at STOP, of shape (~ NUM)."
  (let ((step (/ (- stop start) num)))
    (lazy #'+
     (lazy #'*
      (lazy-index-components (~ num))
      (full (~ num) step))
     (full (~ num) start))))

(defun linspace (start stop &key (num 50) (endpoint 't)  (retstep nil))
  (let* ((step (/ (- stop start) (if endpoint (1- num)  num)))
         (end-ind (if endpoint  num (1- num)))
         (space (lazy #'+
                 (lazy #'*
                  (lazy-index-components (~ end-ind))
                  (full (~ end-ind) step))
                 (full (~ end-ind) start))))
    (if retstep
        `( ,space . ,step)
        space)))

(defun logspace (start stop &key (num 50) (endpoint 't) (base 10.0))
  (lazy (lambda (x) (expt base x))
   (linspace start stop :num num :endpoint endpoint :retstep nil)))



;;; Array manipulation
;;;; Helper functions

(defun negative-shift-down-factor-p (i times axis)
  (and (= i axis) (< times 0)))
(defun positive-shift-down-factor-p (i times axis)
  (and (= i axis) (>= times 0)))




;;;; Exported functions

;; REVIEW: maybe just don't? `lazy-stack' seems perfectly capable here
(defun lazy-concat (arrs &optional (axis 0))
  "Concatenate arrays in list ARRS along AXIS.

Arrays in ARRS must have equal dimensions in all places
except AXIS. Array are harmonized before concatenation.

Numpy has lots of redundant ways to stack/concatenate
arrays, but only one generic one composing `petalisp:lazy-stack' is needed."
  (with-lazy-arrays (arrs)
    (apply #'lazy-stack axis (lazy-harmonize-list-of-arrays arrs))))




;; REVIEW; why not use `peeling-reshaper' here?
;;         for `peeling-reshaper' it is an error to peel
;;         off more layers than the array has, but it might
;;         be fine to allow for the user to shift beyond
;;         the number of array entries
(defun shift-entries (arr &key (axis 0) (times 1) (fill-value 0))
  "Shift elements on axis AXIS of unstrided array away from zero index position;
fill vacated spots with FILL-VALUE and return an unstrided
array of the same shape. TIMES may be negative.

The array MUST be unstrided, i.e. of shape (~ N ~ M ~ K ...) where there
is only one integer after each `~'.

A positive TIMES value would shift away from the zero index, and a negative value
would shift towards the zero index."
  (assert (not (= 0 times)))
  (with-lazy-arrays (arr)
    (let* ((dimensions (shape arr))
           (axis-len (nth axis (shape arr)))
           (fill-shape (apply #'~*
                              (loop for num
                                      in dimensions
                                    for i in (ax:iota (ndim arr))
                                    collect (cond
                                              ((positive-shift-down-factor-p i times  axis)
                                               (range 0 times))
                                              ((negative-shift-down-factor-p i times axis)
                                               (range (+ axis-len times) axis-len))
                                              (t (range num))))))

           (offsets (make-array (ndim arr) :initial-contents (loop repeat (ndim arr) collect 0 ))))
      (setf (aref offsets axis) times)

      (lazy-reshape
       (lazy-fuse-and-harmonize (full  fill-shape fill-value)
                                (lazy-reshape arr (make-transformation
                                                   :offsets offsets )))
       (lazy-array-shape arr)))))



(defun cyclically-shift (arr num axis)
  "Cyclically shift elements of array ARR along AXIS, NUM times."
  (with-lazy-arrays (arr)
    (when (>= axis (length (lazy-array-dimensions arr)))
      (error "Invalid axis ~S for shape ~S." axis (lazy-array-dimensions arr)))
    (let* ((dims (lazy-array-dimensions arr))
           (axis-len (nth axis dims))
           (split-point (mod (- axis-len num) axis-len)))
      (lazy-stack axis (lazy-slices arr (range split-point axis-len) axis)
                  (lazy-slices arr (range 0 split-point) axis)))))


(defun permute-dims (arr &optional (perm nil))
  "Permute axes of ARR according to permutation PERM.

PERM should be a list of indexes from 0..(NDIM ARR) reflecting
the desired final position of the axis number.
If PERM is not specified, reverse axes.

e.g. (permute-dims (ones (~ 2 ~ 3 ~ 4)) '(0 2 1)) will
swap the last two axes."
  (with-lazy-arrays (arr)
    (unless perm
      (setf perm (reverse (ax:iota (ndim arr)))))
    (lazy-reshape arr (petalisp.api:make-transformation
                       :output-mask (make-array
                                     (length perm)
                                     :initial-contents perm)))))



(defun transpose (arr)
  "Return transpose of array.

On 2D arrays, do the matrix transpose.
On ND arrays, reverse axis order."
  (with-lazy-arrays (arr)
    (permute-dims arr)))

(defun matrix-transpose (arr)
  "Matrix transpose a matrix or stack of matrices.

Peforms `swap-axes' on the last two dimensions of ND array, which
should form KxL matrices (i.e shape `(... ~ K ~ L)')"
  (with-lazy-arrays (arr)
    (assert (>= (ndim arr) 2) nil "Cannot matrix transpose an array of rank less than 2.")
    (let ((len (ndim arr))
          (output-mask (make-array (ndim arr)
                                   :initial-contents
                                   (ax:iota (ndim arr)))))
      (setf (aref output-mask (- len 2))  (- len 1))
      (setf (aref output-mask (- len 1)) (- len 2))
      (lazy-reshape arr  (make-transformation
                          :output-mask output-mask)))))

(defun swap-axes (arr axis1 axis2)
  "Interchange axes of ARR.

E.g. on a 2D array, this corresponds to the matrix transpose."
  (with-lazy-arrays (arr)
    (let ((output-mask (make-array (ndim arr)
                                   :initial-contents
                                   (ax:iota (ndim arr)))))
      (setf (aref output-mask axis1) axis2)
      (setf (aref output-mask axis2) axis1)
      (lazy-reshape arr  (make-transformation
                          :output-mask output-mask )))))



(defun flip (arr &key (axis 'all))
  "Flip elements of ARR along AXIS.

AXIS can be an axis number or list of such.
If AXIS is `all', flip elements along all axes."
  (with-lazy-arrays (arr)
    (let* ((axes-to-flip (cond
                           ((integerp axis) `(,axis))
                           ((equal axis 'all)
                            (ax:iota (ndim arr)))))
           (scalings (make-array (ndim arr)
                                 :initial-contents
                                 (loop for axis
                                         below (ndim arr)
                                       collect
                                       (if (member axis axes-to-flip)
                                           -1
                                           1)))))

      (lazy-reshape arr (make-transformation :scalings scalings)
       (lazy-array-shape arr)))))
