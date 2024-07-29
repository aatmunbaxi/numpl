(ql:quickload :petalisp)

(defpackage :numpl
  (:use :common-lisp :petalisp :petalisp.core
   :alexandria)
  (:local-nicknames (nil :petalisp)
                    (:ax :alexandria)))

(in-package :numpl)

;; TODO add type checks to array functions
;;; Array creation routines
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


(defun random-array (shape &optional (limit 256) (element-type nil))
  "Return lazy-array with SHAPE full of random values
with limit LIMIT."
  (let* ((element-type (unless element-type (type-of limit)))
         (array (make-array (shape-dimensions shape) :element-type element-type)))
    (loop for index below (array-total-size array) do
      (setf (row-major-aref array index)
            (random limit)))
    (lazy-array array)))



(defun negative-shift-down-factor-p (i times axis)
  (and (= i axis) (< times 0)))
(defun positive-shift-down-factor-p (i times axis)
  (and (= i axis) (>= times 0)))

;; TODO; Certainly one of the functions of all time...
;; REVIEW; why not use `peeling-reshaper' here?
(defun shift-entries (arr &key (axis 0) (times 1) (fill-value 0))
  "Shift elements on axis AXIS of unstrided array away from zero index position;
fill vacated spots with FILL-VALUE and return an unstrided
array of the same shape. TIMES may be negative.

The array MUST be unstrided, i.e. of shape (~ N ~ M ~ K ...) where there
is only one integer after each `~'.

A positive TIMES value would shift away from the zero index, and a negative value
would shift towards the zero index."
  (with-lazy-arrays (arr)
    (let* ((dimensions (shape arr))
           (axis-len (nth axis (shape arr)))
           (fill-shape   (loop for num in dimensions
                               for i in (iota (ndim arr))
                               append (cond
                                        ((positive-shift-down-factor-p i times  axis)
                                         `(~ 0 ,times))

                                        ((negative-shift-down-factor-p i times axis)
                                         `(~ ,(+ axis-len times) ,axis-len))
                                        (t `(~ ,num)))))
           (axis-labels  (gen-axis-labels arr)))
      (setf (nth axis axis-labels) `( ,(car (nth axis axis-labels)) .
                                      (+ ,times ,(car (nth axis axis-labels)))))
      (lazy-reshape
       (lazy-fuse-and-harmonize (full (eval fill-shape) fill-value)
                                (lazy-reshape arr (eval `(transform ,@(mapcar #'car axis-labels)
                                                                    to
                                                                    ,@(mapcar #'cdr axis-labels)))))
       (lazy-array-shape arr)))))




;; TODO write eye function
(defun eye (N &optional (M nil) (k 0))
  "Return 2D lazy-array with ones on diagnoal, zeros elsewhere."
  (unless M
    (setf M N))
  (let ((max-dim (max N M)))))


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


;; REVIEW: maybe just don't? `lazy-stack' seems perfectly capable here
(defun lazy-concat (arrs &optional (axis 0))
  "Concatenate arrays in list ARRS along AXIS.

Arrays in ARRS must have equal dimensions in all places
except AXIS. Array are harmonized before concatenation.

Numpy has lots of redundant ways to stack/concatenate
arrays, but only one generic one composing `petalisp:lazy-stack' is needed."
  (with-lazy-arrays ('arrs)
    (apply #'lazy-stack axis (lazy-harmonize-list-of-arrays arrs))))




;;; Array manipulation
;;;; Helper functions
(defun gen-permutation (perm)
  "Generate mapping from permutation to be passed to `gen-axis-labels'.
PERM is a list of indexes equinumerous to the number of axes
to permute, where each position in PERM reflects the final position of
the axis."
  (lambda (x) (nth x perm)))

(defun gen-axis-labels (arr &optional (trans 'identity) (prefix "a"))
  "Generate alist of axis labels corresponding to a trans
of axis numberings to be passed to `petalisp.transform'

trans should be a function taking in an integer from 0..N
where N is the rank of ARR which is bijective modulo N.

e.g. to shift the axis labels cyclically to the right by K places,
use `(gen-axis-labels ARR (lambda (i) (+ i k))))'"
  (with-lazy-arrays (arr)
    (let ((num-axes (length (lazy-array-dimensions arr))))
      (loop for index
              below num-axes
            collect
            `(,(intern (concatenate 'string prefix (string
                                                    (format nil "~A" index))))
              .
              ,(intern (concatenate 'string prefix (string
                                                    (format nil "~A" (mod (funcall trans index) num-axes))))))))))

;;;; Exported functions
(defun permute-dims (arr &optional (perm nil))
  "Permute axes of ARR according to permutation PERM.

PERM should be a list of indexes from 0..(NDIM ARR) reflecting
the desired final position of the axis number.
If PERM is not specified, reverse axes.

e.g. (permute-dims (ones (~ 2 ~ 3 ~ 4)) '(0 2 1)) will
swap the last two axes."
  (with-lazy-arrays (arr)
    (unless perm
      (setf perm  (reverse (ax:iota (ndim arr)))))
    (let* ((mapping-alist (gen-axis-labels arr (gen-permutation perm))))
      (lazy-reshape arr  (eval `(transform ,@(mapcar #'car mapping-alist)
                                           to
                                           ,@(mapcar #'cdr mapping-alist)))))))

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
    (flet ((swap-last-two (x) (cond ((= x (1- (ndim arr)))
                                     (1- x))
                                    ((= x (- 2 (ndim arr)))
                                     (1+ x))
                                    (t x))))
      (let ((perm  (gen-axis-labels arr  #'swap-last-two )))
        (lazy-reshape arr  (eval `(transform ,@(mapcar #'car perm) to ,@(mapcar #'cdr perm))))))))

(defun swap-axes (arr axis1 axis2)
  "Interchange axes of ARR.

E.g. on a 2D array, this corresponds to the matrix transpose."
  (with-lazy-arrays (arr)
    (let* ((num-axes (length (lazy-array-dimensions arr)))
           (original-indexes (gen-axis-labels arr))
           (after-indexes (loop for index
                                  below num-axes
                                collect
                                ;; HACK maybe reduce this to one list to avoid iterating twice?
                                (cond ((= index axis1)
                                       (intern (concatenate 'string "a" (format nil "~A" axis2))))
                                      ((= index axis2)
                                       (intern (concatenate 'string "a"  (format nil "~A" axis1))))
                                      (t (intern (concatenate 'string "a"  (format nil "~A" index))))))))
      (lazy-reshape arr  (eval `(transform ,@original-indexes to ,@after-indexes))))))


(defun flip (arr &key (axis 'all))
  "Flip elements of ARR along AXIS.

AXIS can be an axis number or list of such.
If AXIS is `all', flip elements along all axes."
  (with-lazy-arrays (arr)
    (let* ((axes-to-flip (cond
                           ((integerp axis) `(,axis))
                           ((equal axis 'all)
                            (ax:iota (ndim arr)))))
           (axis-labels (gen-axis-labels arr)))
      (loop for pair in axis-labels
            for i in (iota (ndim arr)) do
              (if (member i axes-to-flip)
                  (setf (cdr pair) `(- ,(cdr pair)))))
      (lazy-reshape arr  (eval `(transform ,@(mapcar #'car axis-labels) to ,@(mapcar #'cdr axis-labels)))))))



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
