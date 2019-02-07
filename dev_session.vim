let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
cd ~/Desktop/Science/postdoc_MDC/projects/nanopore
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +1 Snakefile
badd +0 scripts/00_SM/func_defs.py
badd +0 scripts/00_SM/rules_wholefastq.py
badd +0 scripts/00_SM/rules_chunks.py
argglobal
silent! argdel *
$argadd Snakefile
$argadd scripts/00_SM/func_defs.py
set stal=2
tabnew
tabnext -1
edit Snakefile
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
wincmd =
argglobal
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
6,8fold
11,22fold
24,36fold
41,49fold
54,68fold
69,78fold
39,82fold
84,100fold
102,111fold
118,120fold
125,146fold
152,168fold
174,183fold
188,203fold
206,241fold
39
normal! zo
39
normal! zc
let s:l = 118 - ((117 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
118
normal! 0
wincmd w
argglobal
2argu
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
1,4fold
8,9fold
11,17fold
19,21fold
23,32fold
34,44fold
46,57fold
let s:l = 1 - ((0 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
wincmd =
tabnext
edit scripts/00_SM/rules_wholefastq.py
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
wincmd =
argglobal
1argu
if bufexists('scripts/00_SM/rules_wholefastq.py') | buffer scripts/00_SM/rules_wholefastq.py | else | edit scripts/00_SM/rules_wholefastq.py | endif
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
2,13fold
17,30fold
34,45fold
49,63fold
67,81fold
let s:l = 1 - ((0 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
1argu
if bufexists('scripts/00_SM/rules_chunks.py') | buffer scripts/00_SM/rules_chunks.py | else | edit scripts/00_SM/rules_chunks.py | endif
setlocal fdm=manual
setlocal fde=0
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
silent! normal! zE
8,28fold
35,45fold
51,64fold
70,81fold
87,100fold
106,119fold
let s:l = 3 - ((2 * winheight(0) + 22) / 45)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
3
normal! 0
wincmd w
wincmd =
tabnext 1
set stal=1
if exists('s:wipebuf') && len(win_findbuf(s:wipebuf)) == 0
  silent exe 'bwipe ' . s:wipebuf
endif
unlet! s:wipebuf
set winheight=1 winwidth=20 shortmess=filnxtToO
set winminheight=1 winminwidth=1
let s:sx = expand("<sfile>:p:r")."x.vim"
if file_readable(s:sx)
  exe "source " . fnameescape(s:sx)
endif
let &so = s:so_save | let &siso = s:siso_save
doautoall SessionLoadPost
unlet SessionLoad
" vim: set ft=vim :
