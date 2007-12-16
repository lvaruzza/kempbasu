(in-package :kempbasu)


;;;
;;; Function defining Macros 
;;;

(defmacro extract-doc-string! (body)
  `(let ((first-line (first ,body)))
     (if (and (> (length ,body) 1)
	      (stringp (first ,body)))
	 (progn
	   (setf ,body (rest ,body))
	   first-line)
	 "")))

(defun def-IO-fun (name var &key args direction open-options body)
  "Returt the code for a function with NAME of a VAR stream or path, and optinonal extra-args, which execute body.
If the argument is a steam direction and open-options are used in function with-open-file" 
  (let ((fun (gensym))
	(stream (gensym))
	(doc-string (extract-doc-string! body)))    
    `(defun ,name (,var ,@args)
       ,doc-string
       (flet ((,fun (,var)
		,@body))
	 (if (or (stringp ,var) (pathnamep ,var))
	     (with-open-file (,stream ,var :direction ,direction ,@open-options)
	       (,fun ,stream))
	     (,fun ,var))))))


(defmacro def-input-fun (name (var &rest args) &body body)
  "def-IO-fun for input file/steam"
  (def-IO-fun name var :args args :direction :input :body body))

(defmacro def-ouput-fun (name (var &rest args) &body body)
  "def-IO-fun for output file/steam"
  (def-IO-fun name var :args args :direction :output
	      :open-options '(:if-exists :supersede)
	      :body body))

;;;
;;; Line oriented Processing
;;;

(defmacro each-line (var stream &body body)
  "Execute body for each line of stream"
  `(loop for ,var = (read-line ,stream nil nil)
      while ,var
      do ,body))

(defun collect-lines (stream f)
  "Collect (f line) for the lines  of stream"
  (loop for line = (read-line stream nil nil)
     while line
     collect (funcall f line)))

(defmacro each-line (var stream &body body)
  "Execute body for each line of stream"
  `(loop for ,var = (read-line ,stream nil nil)
      while ,var
      do ,body))
