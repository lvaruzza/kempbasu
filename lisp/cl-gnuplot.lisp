;Ryans gnuplot interface library
; http://www.nonperiodic.net/blog/2005/07/22/more-gnuplotting-lispers/
; http://www.cliki.net/Common%20Lisp%20and%20gnuplot

; Modifications added for port to SBCL on OSX
; By Joseph Lee Nunn III
; 10/05/2006

(defpackage :cl-gnuplot
  (:use :common-lisp)
  (:export #:test-plot
	   #:test-animation))

(in-package :cl-gnuplot)

(defun string-join (delimiter strings)
  (format nil (format nil "~~{~~a~~^~a~~}" delimiter) strings))

;;; SERIES-2D CLASS
(defclass series-2d ()
  ((data :initarg :data))
  (:documentation
   "A class for managing data to be plotted with gnuplot."))

(defun make-series-2d (data)
  (make-instance 'series-2d :data data))

;;; PLOT CLASS
(defclass plot ()
  ((series    :initform ()
              :allocation :instance)
   (seriesf   :initform (make-hash-table))
   
   (title     :initform nil)
   (xlabel    :initform nil)
   (x2label   :initform nil)
   (ylabel    :initform nil)
   (y2label   :initform nil)
    
   (process   :initform 
	      #+sbcl (sb-ext:run-program;extensions:run-program
				"gnuplot"        ; program to run
				'("-display" ":0.0")  ; need this non nil for osx to work
				:wait   nil      ; wait for child
				:pty    nil      ; no terminal
				:search t
				:input  :stream  ; need an input stream
				:error  :output  ; print to stdout for debugging
				:output t)       ; print to stdout for debugging 
	      #+cmu (extensions:run-program
				"gnuplot"        ; program to run
				nil              ; no arguments
				:wait   nil      ; wait for child
				:pty    nil      ; no terminal
				:input  :stream  ; need an input stream
				:error  :output  ; print to stdout for debugging
				:output t)       ; print to stdout for debugging
              :allocation :instance)   
   (stream    :initform nil))
  (:documentation
   "A class for plotting data with gnuplot."))

(defmethod plot-title ((p plot) &optional value)
  (unless (null value)
    (setf (slot-value p 'title) value)
    (format (slot-value p 'stream) "set title \"~a\"~%" value))
  (slot-value p 'title))
(defmethod plot-xlabel ((p plot) &optional value)
  (unless (null value)
    (setf (slot-value p 'xlabel) value)
    (format (slot-value p 'stream) "set xlabel \"~a\"~%" value))
  (slot-value p 'xlabel))
(defmethod plot-x2label ((p plot) &optional value)
  (unless (null value)
    (setf (slot-value p 'x2label) value)
    (format (slot-value p 'stream) "set x2label \"~a\"~%" value))
  (slot-value p 'x2label))   
(defmethod plot-ylabel ((p plot) &optional value)
  (unless (null value)
    (setf (slot-value p 'ylabel) value)
    (format (slot-value p 'stream) "set ylabel \"~a\"~%" value))
  (slot-value p 'ylabel))
(defmethod plot-y2label ((p plot) &optional value)
  (unless (null value)
    (setf (slot-value p 'y2label) value)
    (format (slot-value p 'stream) "set y2label \"~a\"~%" value))
  (slot-value p 'y2label))

(defun make-plot ()
  (make-instance 'plot))

(defmethod initialize-instance :after ((plt plot) &key)
  (setf (slot-value plt 'stream)
        #+sbcl (sb-ext:process-input (slot-value plt 'process))
	#+cmu (extensions:process-input (slot-value plt 'process))))

(defmethod plot-close ((plt plot))
  #+sbcl (sb-ext:process-close (slot-value plt 'process))
  #+cmu (extensions:process-close (slot-value plt 'process)))

(defmethod plot-cmd ((plt plot) command)
  (format (slot-value plt 'stream) "~a~%" command))

(defmethod plot-draw ((plt plot))
  (let ((plot-string ()))

   ; Loop through the series, sending format information
    (dolist (subplot (slot-value plt 'series))
      (let ((subplot-string '("'-'"))
            (sf (gethash subplot (slot-value plt 'seriesf))))
        (unless (null (series-format-axes sf))
          (push "axes" subplot-string)
          (push (series-format-axes sf) subplot-string))
        (unless (null (series-format-title sf))
          (push "title" subplot-string)
          (push (format nil "\"~a\"" (series-format-title sf))
                subplot-string))
        (unless (null (series-format-style sf))
          (push "with" subplot-string)
          (push (series-format-style sf) subplot-string))
        (unless (null (series-format-smooth sf))
          (push "smooth" subplot-string)
          (push (series-format-smooth sf) subplot-string))
        (push (string-join " " (reverse subplot-string)) plot-string)))

    ; Join our list of subplots and print them to the stream
    (format (slot-value plt 'stream) "plot ~a~%"
            (string-join ", " (reverse plot-string)))

    ; Loop through the series, sending the data
    (dolist (subplot (slot-value plt 'series))
      (dolist (line (slot-value subplot 'data))
        (format (slot-value plt 'stream) "~{~a~^ ~}~%" line))
      (format (slot-value plt 'stream) "e~%"))
    
    ; Make sure the plot string is written to the stream
    (force-output (slot-value plt 'stream))))

(defgeneric plot-add-series (plot series-2d &key))

(defmethod plot-add-series
    ((plt plot) (series series-2d) &key style title smooth axes)
  (push series (slot-value plt 'series))
  (let ((sf (make-series-format)))
    (setf (gethash series (slot-value plt 'seriesf)) sf)
    (unless (null style)  (setf (series-format-style  sf) style))
    (unless (null title)  (setf (series-format-title  sf) title))
    (unless (null axes)   (setf (series-format-axes   sf) axes))
    (unless (null smooth) (setf (series-format-smooth sf) smooth))))

(defmethod plot-remove-series
    ((p plot) (series series-2d))
  (remhash series (slot-value p 'seriesf))
  (setf (slot-value p 'series) (remove series (slot-value p 'series))))

(defstruct series-format
  style
  title
  smooth
  axes)

(defun test-plot ()
  (let ((series1
         (make-series-2d '((5 1) (6 -2) (7 4) (8 -8) (9 16) (12 -32))))
        (series2
         (make-series-2d '((5) (6) (8) (2) (3) (8))))
        (series3
         (make-series-2d
          '((1 2) (2 3) (3 5) (4 7) (5 11) (6 13) (7 17))))
        (plot (make-plot)))
    (plot-add-series plot series1
                     :style  "linespoints"
                     :title  "Foo"
                     :axes   "x2y1")
    (plot-add-series plot series2
                     :style  "lines"
                     :title  "Bar"
                     :smooth "bezier"
                     :axes   "x2y1")
    (plot-add-series plot series3
                     :style  "points"
                     :title  "Primes"
                     :axes   "x1y2")

    (plot-title plot "Test Title")
    (plot-xlabel plot "X1 Label")
    (plot-x2label plot "X2 Label")
    (plot-ylabel plot "Y1 Label")
    (plot-y2label plot "Y2 Label")

    (plot-cmd plot "set logscale y")
    
    (plot-draw plot)

    (break "Write postscript file?")
     
    (plot-close plot)))

(defun test-animation ()
  (let ((plot (make-plot)))

    ; Start each frame at a different spot
    (dotimes (count 1000)
      (let ((start (* 0.1 count))
            (data ()))

        ; Generate each data point
        (dotimes (len 100)
          (push (list (sin (+ start (* len 0.1)))) data))

        ; Make the series
        (let ((series (make-series-2d data)))
          (plot-add-series plot series :style "lines")
          (plot-draw plot)
          (plot-remove-series plot series))
        (sleep 0.05)))
    (break "Close gnuplot?")
    (plot-close plot)))
