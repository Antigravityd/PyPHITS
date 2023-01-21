;;; Phits-mode.el --- Summary
;;; Commentary:
;;

;;; Code:

(require 'font-lock)
(require 'align)


(defvar phits-mode-syntax-table
  (let ((st (make-syntax-table))) ;; borrowed from fortran-mode.el syntax table
    (modify-syntax-entry ?$  "<"  st)
    (modify-syntax-entry ?\; "."  st)
    (modify-syntax-entry ?\r " "  st)
    (modify-syntax-entry ?+  "."  st)
    (modify-syntax-entry ?-  "."  st)
    (modify-syntax-entry ?=  "."  st)
    (modify-syntax-entry ?*  "."  st)
    (modify-syntax-entry ?/  "."  st)
    (modify-syntax-entry ?\' "\"" st)
    (modify-syntax-entry ?\" "\"" st)
    (modify-syntax-entry ?\\ "\\" st)
    (modify-syntax-entry ?.  "_"  st)
    (modify-syntax-entry ?_  "_"  st)
    (modify-syntax-entry ?\! "<"  st)
    (modify-syntax-entry ?\# "<"  st) ;; allowing "#" as a comment mark will break on some [Cell] sections
    (modify-syntax-entry ?\n ">"  st)
    (modify-syntax-entry ?%  "<"  st)
    st)
  "Syntax table used for .inp files.")

;; Line continuation isn't handled at all

(setq phits-comment-regexp "^ \\{,4\\}c.*") ;; for the archaic fixed-format comments
(setq phits-section-regexp "^ \\{,4\\}\\[.*\\]")
(setq phits-parameter-regexp "\\s-*\\([[:alnum:]<>-]+\\)\\((.*)\\|\\[.*\\]\\)?\\s-*=")
(setq phits-label-regexp "^\\s-*\\(\\w*\\):")

(setq phits-function-regexp
      (regexp-opt
       '("float" "int" "abs""exp" "log" "log10" "max" "min" "mod" "nint" "sign" "sqrt" "acos" "asin" "atan"
	 "atan2" "cos" "cosh" "sin" "sinh" "tan" "tanh" "pi")
       'words))


(setq phits-particle-regexp
      (concat "\\(\\<\\(\\([0-9]\\{,3\\}"
	      (regexp-opt '("H" "He" "Li" "Be" "B" "C" "N" "O" "F" "Ne" "Na" "Mg" "Al" "Si" "P" "S" "Cl" "Ar"
			    "K" "Ca" "Sc" "Ti" "V" "Cr" "Mn" "Fe" "Co" "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr"
			    "Rb" "Sr" "Y" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh" "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "I" "Xe"
			    "Cs" "Ba" "La" "Ce" "Pr" "Nd" "Pm" "Sm" "Eu" "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu"
			    "Hf" "Ta" "W" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl" "Pb" "Bi" "Po" "At" "Rn" "Fr" "Ra" "Ac"
 			    "Th" "Pa" "U" "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es" "Fm" "Md" "No" "Lr" "Rf" "Db" "Sg" "Bh"
			    "Hs" "Mt" "Ds" "Rg" "Cn" "Nh" "Fl" "Mc" "Lv" "Ts" "Og")
			  t)
	      "\\)\\|\\("
	      (regexp-opt '("H" "He" "Li" "Be" "B" "C" "N" "O" "F" "Ne" "Na" "Mg" "Al" "Si" "P" "S" "Cl" "Ar"
			    "K" "Ca" "Sc" "Ti" "V" "Cr" "Mn" "Fe" "Co" "Ni" "Cu" "Zn" "Ga" "Ge" "As" "Se" "Br" "Kr"
			    "Rb" "Sr" "Y" "Zr" "Nb" "Mo" "Tc" "Ru" "Rh" "Pd" "Ag" "Cd" "In" "Sn" "Sb" "Te" "I" "Xe"
			    "Cs" "Ba" "La" "Ce" "Pr" "Nd" "Pm" "Sm" "Eu" "Gd" "Tb" "Dy" "Ho" "Er" "Tm" "Yb" "Lu"
			    "Hf" "Ta" "W" "Re" "Os" "Ir" "Pt" "Au" "Hg" "Tl" "Pb" "Bi" "Po" "At" "Rn" "Fr" "Ra" "Ac"
 			    "Th" "Pa" "U" "Np" "Pu" "Am" "Cm" "Bk" "Cf" "Es" "Fm" "Md" "No" "Lr" "Rf" "Db" "Sg" "Bh"
			    "Hs" "Mt" "Ds" "Rg" "Cn" "Nh" "Fl" "Mc" "Lv" "Ts" "Og")
			  t)
	      "\\(-[0-9]\\{,3\\}\\)\\)\\|" ;; collisions common for "all," "other," and "pi"
	      (regexp-opt '("all" "proton" "neutron" "pion0"  "kaon0" "other" "electron" "positron"  "photon"
			    "gamma" "deuteron" "triton" "3he" "alpha" "nucleus" "pi")
			  t)
	      "\\)\\>\\|pion\\+\\|pion-\\|kaon\\+\\|kaon-\\|muon\\+\\|muon-\\)"))


(setq phits-special-regexp
      (concat "\\<\\(mat\\[.*\\]\\|mt?[0-9]+\\|tr[0-9]+\\|"
	      (regexp-opt '("like" "but" "p" "px" "py" "pz" "so" "s" "sx" "sy" "sz" "c/x" "c/y" "c/z" "cx" "cy" "cz"
			    "k/x" "k/y" "k/z" "kx" "ky" "kz" "sq" "gq" "tx" "ty" "tz" "box" "rpp" "sph" "rcc" "rhp"
			    "hex" "rec" "trc" "ell" "wed")
			  t)
	      "\\)\\>"))

;; custom matchers required for selectively case-insensitive fontification
(defun phits-match-archaic-comment (limit)
  (let ((case-fold-search t))
    (re-search-forward phits-comment-regexp limit t)))

(defun phits-match-function (limit)
  (let ((case-fold-search t))
    (re-search-forward phits-function-regexp limit t)))

(defun phits-match-special (limit)
  (let ((case-fold-search t))
    (re-search-forward phits-special-regexp limit t)))

(defun phits-line-matching (regexp)
    (save-excursion
       (beginning-of-line)
       (if (re-search-forward regexp (line-end-position) t)
	   t)))

(defun phits-match-nongrid ()
  (not (or
       (phits-line-matching phits-section-regexp)
       (phits-line-matching phits-parameter-regexp)
       (phits-line-matching phits-label-regexp)
       (phits-line-matching "^\\(\\s-*\\)\\w+\\[.*\\]"))))

(defvar phits-archaic-comment-font-lock ;; disabled by default; add this to phits-font-lock if you don't mind it turning your lines starting with Co in [Material] the wrong color
  (cons 'phits-match-archaic-comment  font-lock-comment-face))

(defvar phits-section-font-lock
  (cons phits-section-regexp font-lock-warning-face))

(defvar phits-parameter-font-lock
  (list phits-parameter-regexp 1 font-lock-variable-name-face))

(defvar phits-label-font-lock
  (list phits-label-regexp 1 font-lock-keyword-face))

(defvar phits-particle-font-lock
  (list phits-particle-regexp 1 font-lock-function-name-face))

(defvar phits-function-font-lock
  (cons 'phits-match-function  font-lock-builtin-face))

;; TODO: solve collision of p, s, and u with the isotopes above
(defvar phits-special-font-lock
  (cons 'phits-match-special font-lock-type-face))


(setq phits-font-lock
      (list phits-section-font-lock
	    phits-parameter-font-lock
	    phits-label-font-lock
	    phits-particle-font-lock
	    phits-function-font-lock
	    phits-special-font-lock))

(defun phits-indent-line () ;; The only thing this fails on is gridlike sections inside a block of parameters,
  (align-region nil
		nil
		'group
		'((phits-parameter
		   (regexp   . "^\\(\\s-*[[:alnum:]<>()-]+\\(?:(.*)\\|\\[.*\\]\\)?\\)\\(\\s-*\\)=\\(\\s-*\\)")
		   (group    . (1 2 3))
		   (justify  . t)
		   (separate . "^\\s-*\\[.*\\]"))
		  (phits-label
		   (regexp   . "^\\(\\s-*\\)\\w+\\(\\s-*\\):\\(\\s-*\\)")
		   (group    . (1 3))
		   (spacing  . (0 1)))
		  (phits-material
		   (regexp   . "^\\(\\s-*\\)\\w+\\[.*\\]")
		   (spacing  . 0))
		  (phits-grid
		   (run-if   . phits-match-nongrid)
		   (regexp   . "\\(\\s-*\\)[[:alnum:]()*/+.-]+")
		   (group    . 1)
		   (repeat   . t)
		   (valid    . phits-match-nongrid))
		  (phits-comment
		   (regexp   . "\\S-\\(\\s-*\\)[#!$%]")))
		'((phits-exclude-after-comment
		   (regexp   . "\\([#!$%].*\\)")))))


(defvar phits-command "phits"
  "The shell command through which PHITS is to be run.  Should be everything but the filename (which will be placed at the end).")

(defvar phits-eps-viewer "gv"
  "The shell command for launching the viewer for any .eps files produced by the simulation.")

(defvar phits-eps-options "--orientation='landscape' --presentation"
  "String of options passed to the .eps viewer.")

(defun run-phits ()
  "Run phits-command on the current buffer asynchronously, displaying output in a new buffer and placing output files in the same directory as the current buffer's file."
  (interactive)
  (let ((default-directory (file-name-directory buffer-file-name))
	(file (buffer-file-name)))
    (switch-to-buffer-other-window (get-buffer-create "PHITS Output"))
    (end-of-buffer)
    (insert (concat (current-time-string) "\n"))
    (start-process-shell-command "run-phits" "PHITS Output" (concat phits-command " " file))))

(defun async-shell-command-no-window (command)
  "Calls command, but doesn't make a buffer containing its output."
  (interactive (list (read-shell-command "$ ")))
  (start-process-shell-command command nil command))

(defun phits-view-results ()
  "Open summary of tally results in external image viewer."
  (interactive)
  (let* ((default-directory (file-name-directory buffer-file-name))
	 (choices (directory-files default-directory nil "\\.eps"))
	 (file (completing-read "View results of tally: " choices))
	 (async-shell-command-display-buffer nil))
    (async-shell-command-no-window (concat phits-eps-viewer " " phits-eps-options " " file))))

(defun phits-check-geometry ()
  "Copies all geometry definitions to a test .inp file with inctl=8 and a proper tally with gshow set.")

(defun phits-abort-run ()
  "Kills current run-phits process, if any."
  (interactive)
  (if (get-process "run-phits")
      (kill-process "run-phits")))

;;;###autoload
(define-derived-mode phits-mode prog-mode "PHITS Input"
  "Testing mode I whipped up that's inspired in small part by https://github.com/kbat/mc-tools/blob/master/mctools/phits/phits-mode.el"
  :syntax-table phits-mode-syntax-table
  (setq-local font-lock-defaults `(,phits-font-lock nil nil))
  (setq-local indent-line-function #'phits-indent-line)
  (add-hook 'phits-mode-hook
	    (lambda ()
	      (local-set-key (kbd "C-c C-c") #'run-phits)
	      (local-set-key (kbd "C-c C-v") #'phits-view-results)
	      (local-set-key (kbd "C-c C-k") #'phits-abort-run)
	      (local-set-key (kbd "C-c C-i") #'phits-check-geometry))))

(provide 'phits-mode)

;;; phits-mode.el ends here
