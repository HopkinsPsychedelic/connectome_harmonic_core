git add -A
git commit -m 'automated commit'
git push
docker build -t winstonian3/connectome_harmonic --build-arg SSH_KEY="$(cat ~/.ssh/id_rsa)" --build-arg CACHE_DATE="$(date)" .
docker images
