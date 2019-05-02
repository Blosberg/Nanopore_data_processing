let SessionLoad = 1
if &cp | set nocp | endif
let s:so_save = &so | let s:siso_save = &siso | set so=0 siso=0
let v:this_session=expand("<sfile>:p")
silent only
silent tabonly
if expand('%') == '' && !&modified && line('$') <= 1 && getline(1) == ''
  let s:wipebuf = bufnr('%')
endif
set shortmess=aoO
badd +207 nanopiper.py
badd +4 Snakefile
badd +174 scripts/00_SM/func_defs.py
badd +5 config.json
badd +0 scripts/00_SM/rules_chunks.py
argglobal
silent! argdel *
$argadd nanopiper.py
$argadd Snakefile
set stal=2
tabnew
tabnext -1
edit nanopiper.py
set splitbelow splitright
wincmd _ | wincmd |
vsplit
1wincmd h
wincmd w
wincmd _ | wincmd |
split
1wincmd k
wincmd w
wincmd t
set winminheight=0
set winheight=1
set winminwidth=0
set winwidth=1
wincmd =
argglobal
setlocal fdm=expr
setlocal fde=FoldPythonMain()
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
21
normal! zo
32
normal! zo
140
normal! zo
163
normal! zo
172
normal! zo
231
normal! zo
255
normal! zo
let s:l = 258 - ((46 * winheight(0) + 25) / 50)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
258
normal! 066|
wincmd w
argglobal
2argu
if bufexists('scripts/00_SM/func_defs.py') | buffer scripts/00_SM/func_defs.py | else | edit scripts/00_SM/func_defs.py | endif
setlocal fdm=expr
setlocal fde=FoldPythonFuncdefs()
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
76
normal! zo
let s:l = 124 - ((18 * winheight(0) + 12) / 25)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
124
normal! 0
wincmd w
argglobal
2argu
if bufexists('config.json') | buffer config.json | else | edit config.json | endif
setlocal fdm=expr
setlocal fde=FoldSnakemake()
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 2 - ((1 * winheight(0) + 12) / 24)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
2
normal! 02|
wincmd w
3wincmd w
wincmd =
tabnext
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
1argu
if bufexists('Snakefile') | buffer Snakefile | else | edit Snakefile | endif
setlocal fdm=expr
setlocal fde=FoldSnakemake()
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 1 - ((0 * winheight(0) + 25) / 50)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
1
normal! 0
wincmd w
argglobal
1argu
if bufexists('scripts/00_SM/rules_chunks.py') | buffer scripts/00_SM/rules_chunks.py | else | edit scripts/00_SM/rules_chunks.py | endif
setlocal fdm=expr
setlocal fde=FoldSnakemake()
setlocal fmr={{{,}}}
setlocal fdi=#
setlocal fdl=0
setlocal fml=1
setlocal fdn=20
setlocal fen
let s:l = 10 - ((9 * winheight(0) + 25) / 50)
if s:l < 1 | let s:l = 1 | endif
exe s:l
normal! zt
10
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
