export GUIX_PROFILE=/home/bosberg/projects/nanopore/guix/.guix-profile

export GUIX_LOCPATH=${GUIX_PROFILE}/lib/locale
export GIT_SSL_CAINFO=$GUIX_PROFILE/etc/ssl/certs/ca-certificates.crt
source ${GUIX_PROFILE}/etc/profile   # now the environmental variables will be augmented (from the front)

unset PYTHONPATH

