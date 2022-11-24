filetype plugin indent on
syntax on

syntax reset
set nowrap
"endif

syntax match aaA "[aA]" 
syntax match aaR "[rR]" 
syntax match aaN "[nN]" 
syntax match aaD "[dD]" 
syntax match aaC "[cC]" 
syntax match aaQ "[qQ]" 
syntax match aaE "[eE]" 
syntax match aaG "[gG]" 
syntax match aaH "[hH]" 
syntax match aaI "[iI]" 
syntax match aaL "[lL]" 
syntax match aaK "[kK]" 
syntax match aaM "[mM]" 
syntax match aaF "[fF]" 
syntax match aaP "[pP]" 
syntax match aaS "[sS]" 
syntax match aaT "[tT]" 
syntax match aaW "[wW]" 
syntax match aaY "[yY]" 
syntax match aaV "[vV]" 
syntax match aaB "[bB]" 
syntax match aaX "[xX]" 
syntax match aaZ "[zZ]" 
syntax match aaGap "[/\-/\.]" 
syntax match header "^.\{-}\t"
" set Line Numbering color to gray
highlight LineNr ctermfg=grey

" Amino Acid Coloring 
highlight	aaA	cterm=NONE	guibg=NONE	guifg=#0087ff   	ctermbg=NONE	ctermfg=33
highlight	aaR	cterm=NONE	guibg=NONE	guifg=#af0000    	ctermbg=NONE	ctermfg=124
highlight	aaN	cterm=NONE	guibg=NONE	guifg=#00af00   	ctermbg=NONE	ctermfg=34
highlight	aaD	cterm=NONE	guibg=NONE	guifg=#af00af   	ctermbg=NONE	ctermfg=127
highlight	aaC	cterm=NONE	guibg=NONE	guifg=#d75f87  	ctermbg=NONE	ctermfg=168
highlight	aaQ	cterm=NONE	guibg=NONE	guifg=#00af00   	ctermbg=NONE	ctermfg=34
highlight	aaE	cterm=NONE	guibg=NONE	guifg=#af00af   	ctermbg=NONE	ctermfg=127
highlight	aaG	cterm=NONE	guibg=NONE	guifg=#d75f00   	ctermbg=NONE	ctermfg=166
highlight	aaH	cterm=NONE	guibg=NONE	guifg=#00afaf   	ctermbg=NONE	ctermfg=37
highlight	aaI	cterm=NONE	guibg=NONE	guifg=#0087ff   	ctermbg=NONE	ctermfg=33
highlight	aaL	cterm=NONE	guibg=NONE	guifg=#0087ff   	ctermbg=NONE	ctermfg=33
highlight	aaK	cterm=NONE	guibg=NONE	guifg=#af0000    	ctermbg=NONE	ctermfg=124
highlight	aaM	cterm=NONE	guibg=NONE	guifg=#0087ff   	ctermbg=NONE	ctermfg=33
highlight	aaF	cterm=NONE	guibg=NONE	guifg=#0087ff   	ctermbg=NONE	ctermfg=33
highlight	aaP	cterm=NONE	guibg=NONE	guifg=#d7af00   	ctermbg=NONE	ctermfg=178
highlight	aaS	cterm=NONE	guibg=NONE	guifg=#00af00   	ctermbg=NONE	ctermfg=34
highlight	aaT	cterm=NONE	guibg=NONE	guifg=#00af00   	ctermbg=NONE	ctermfg=34
highlight	aaW	cterm=NONE	guibg=NONE	guifg=#0087ff   	ctermbg=NONE	ctermfg=33
highlight	aaY	cterm=NONE	guibg=NONE	guifg=#00afaf   	ctermbg=NONE	ctermfg=37
highlight	aaV	cterm=NONE	guibg=NONE	guifg=#0087ff   	ctermbg=NONE	ctermfg=33
highlight	aaB	cterm=NONE	guibg=NONE	guifg=NONE 	ctermbg=NONE	ctermfg=NONE
highlight	aaX	cterm=NONE	guibg=NONE	guifg=NONE 	ctermbg=NONE	ctermfg=NONE
highlight	aaZ	cterm=NONE	guibg=NONE	guifg=NONE 	ctermbg=NONE	ctermfg=NONE
highlight	header	cterm=NONE	guibg=NONE	guifg=NONE 	ctermbg=NONE	ctermfg=NONE

set runtimepath=$VIMMOVE
",$VIMRUNTIME


" Vim move keeys
let g:move_key_modifier = 'C'
let g:move_key_modifier_visualmode = 'S'
