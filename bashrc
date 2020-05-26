PS1='\[\033[01;34m\][\W] \[\033[01;32m\]€\[\033[00m\] '
alias lsa="ls -alh"

HISTCONTROL=ignoredups:erasedups
shopt -s histappend
PROMPT_COMMAND="history -n; history -w; history -c; history -r; $PROMPT_COMMAND"

jsgget() {
        f=https://raw.githubusercontent.com/joonatansipola/stuff/master/$1
        wget $f
}
