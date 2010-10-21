;; To load the following function into Emacs : M-x load-file
;; This function only operates on a region. To run this function,
;; first select a region and call M-x extract-header
(defun extract-header ()
  (interactive)
  "Extract the prototypes of the functions and convert them into the documentation format"
  (if (not (and transient-mark-mode mark-active))
    (error "%" "No region active"))
  (let* ((start (region-beginning))
          (end (region-end))
          (line-start (line-number-at-pos start))
          (line-end (line-number-at-pos end)))
    (replace-regexp 
      "extern +\\([a-zA-Z0-9_ \\*]+\\) +\\([a-zA-Z0-9_]+\\) *(\\(.*\\));"
      "\\\\item \\\\describefun{\\1}{\\2}{\\3}"
      nil start end)
    (goto-char (point-min))
    (let ((start (point-at-bol line-start))
           (end (point-at-bol (1+ line-end))))
      (replace-string "*" "\\ptr " nil start end))
    (goto-char (point-min))
    (let ((start (point-at-bol line-start))
           (end (point-at-bol (1+ line-end))))
      (replace-regexp "\\(Pnl[a-zA-Z]*\\)" "\\\\refstruct{\\1}" nil start end)
        )
      )
    )


