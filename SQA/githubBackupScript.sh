#!/bin/bash 
# A simple script to backup an organization's GitHub repositories.
 
GHBU_BACKUP_DIR=${GHBU_BACKUP_DIR-"/pathTo/GithubBackups"}           # where to place the backup files
GHBU_ORG=${GHBU_ORG-"NaluCFD"}                                       # the GitHub organization whose repos will be backed up
                                                                     # (if you're backing up a user's repos instead, this should be your GitHub username)
GHBU_UNAME=${GHBU_UNAME-"spdomin"}                                   # the username of a GitHub account (to use with the GitHub API)
GHBU_GITHOST=${GHBU_GITHOST-"github.com"}                            # the GitHub hostname (see comments)
GHBU_PRUNE_OLD=${GHBU_PRUNE_OLD-true}                                # when `true`, old backups will be deleted
GHBU_PRUNE_AFTER_N_DAYS=${GHBU_PRUNE_AFTER_N_DAYS-7}                 # the min age (in days) of backup files to delete
GHBU_SILENT=${GHBU_SILENT-false}                                     # when `true`, only show error messages 
GHBU_API=${GHBU_API-"https://api.github.com"}                        # base URI for the GitHub API
GHBU_GIT_CLONE_CMD="git clone --mirror git@${GHBU_GITHOST}:"         # base command to use to clone GitHub repos
REPOLIST='Nalu NaluDoc NaluSandBox'                                  # list of desired products under user name (avoids password usage)

TSTAMP=`date "+%Y%m%d-%H%M"`
 
# The function `check` will exit the script if the given command fails.
function check {
  "$@"
  status=$?
  if [ $status -ne 0 ]; then
    echo "ERROR: Encountered error (${status}) while running the following:" >&2
    echo "           $@"  >&2
    echo "       (at line ${BASH_LINENO[0]} of file $0.)"  >&2
    echo "       Aborting." >&2
    exit $status
  fi
}
 
# The function `tgz` will create a gzipped tar archive of the specified file ($1) and then remove the original
function tgz {
   check tar zcf $1.tar.gz $1 && check rm -rf $1
}

$GHBU_SILENT || (echo "" && echo "=== INITIALIZING ===" && echo "")
 
$GHBU_SILENT || echo "Using backup directory $GHBU_BACKUP_DIR"
check mkdir -p $GHBU_BACKUP_DIR

$GHBU_SILENT || (echo "" && echo "=== BACKING UP ===" && echo "")

for REPO in $REPOLIST; do
$GHBU_SILENT || echo "Backing up ${GHBU_ORG}/${REPO}"
check ${GHBU_GIT_CLONE_CMD}${GHBU_ORG}/${REPO}.git ${GHBU_BACKUP_DIR}/${GHBU_ORG}-${REPO}-${TSTAMP}.git && tgz ${GHBU_BACKUP_DIR}/${GHBU_ORG}-${REPO}-${TSTAMP}.git
 
$GHBU_SILENT || echo "Backing up ${GHBU_ORG}/${REPO}.wiki"
${GHBU_GIT_CLONE_CMD}${GHBU_ORG}/${REPO}.wiki.git ${GHBU_BACKUP_DIR}/${GHBU_ORG}-${REPO}.wiki-${TSTAMP}.git 2>/dev/null && tgz ${GHBU_BACKUP_DIR}/${GHBU_ORG}-${REPO}.wiki-${TSTAMP}.git

$GHBU_SILENT || echo "Backing up ${GHBU_ORG}/${REPO}.issues"
${GHBU_GIT_CLONE_CMD}${GHBU_ORG}/${REPO}/issues ${GHBU_BACKUP_DIR}/${GHBU_ORG}-${REPO}.issues-${TSTAMP}.git 2>/dev/null && tgz ${GHBU_BACKUP_DIR}/${GHBU_ORG}-${REPO}.issues-${TSTAMP}.git
done
 
if $GHBU_PRUNE_OLD; then
  $GHBU_SILENT || (echo "" && echo "=== PRUNING ===" && echo "")
  $GHBU_SILENT || echo "Pruning backup files ${GHBU_PRUNE_AFTER_N_DAYS} days old or older."
  $GHBU_SILENT || echo "Found `find $GHBU_BACKUP_DIR -name '*.tar.gz' -mtime +$GHBU_PRUNE_AFTER_N_DAYS | wc -l` files to prune."
  find $GHBU_BACKUP_DIR -name '*.tar.gz' -mtime +$GHBU_PRUNE_AFTER_N_DAYS -exec rm -fv {} > /dev/null \; 
fi
 
$GHBU_SILENT || (echo "" && echo "=== DONE ===" && echo "")
$GHBU_SILENT || (echo "GitHub backup completed." && echo "")
