#' ---
#' title: googledrive test clean
#' date: '`r format(Sys.time())`'
#' output: html_document
#' ---

#' This script aggregates the test-related clean code from all test files.

library(googledrive)
source(testthat::test_path('helper.R'))
whoami <- drive_user()$user
whoami[c('displayName', 'emailAddress')]

## change this to TRUE when you are really ready to do this!
CLEAN <- FALSE

#' ## test-drive_cp.R
me_ <- nm_fun("TEST-drive-cp")
nm_ <- nm_fun("TEST-drive-cp", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("i-am-a-folder"),
    nm_("not-unique-folder"),
    nm_("i-am-a-file")
  ))
}
#' ## test-drive_create.R
me_ <- nm_fun("TEST-drive-create")
nm_ <- nm_fun("TEST-drive-create", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("create-in-me"),
    nm_("create-me-in-root")
  ))
}
#' ## test-drive_download.R
nm_ <- nm_fun("TEST-drive-download", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("DESC"),
    nm_("DESC-doc")
  ))
}
#' ## test-drive_find.R
me_ <- nm_fun("TEST-drive-find")
nm_ <- nm_fun("TEST-drive-find", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("find-me"),
    nm_("this-should-not-exist")
  ))
}
#' ## test-drive_get.R
nm_ <- nm_fun("TEST-drive-get", NULL)
if (CLEAN) {
  files <- drive_find(nm_("thing0[1234]"))
  drive_trash(files)
  parents <- drive_find(nm_("parent0[12]"))
  drive_trash(parents)
  drive_trash(nm_("child_of_2_parents"))
}
#' ## test-drive_ls.R
nm_ <- nm_fun("TEST-drive-ls", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("list-me"),
    nm_("this-should-not-exist")
  ))
}
#' ## test-drive_mkdir.R
me_ <- nm_fun("TEST-drive-mkdir")
nm_ <- nm_fun("TEST-drive-mkdir", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("OMNI-PARENT"),
    nm_("I-live-in-root")
  ))
}
#' ## test-drive_mv.R
me_ <- nm_fun("TEST-drive-mv")
nm_ <- nm_fun("TEST-drive-mv", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("move-files-into-me"),
    nm_("DESC"),
    nm_("DESC-renamed")
  ))
}
#' ## test-drive_publish.R
nm_ <- nm_fun("TEST-drive-publish", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("foo_pdf"),
    nm_("foo_doc"),
    nm_("foo_sheet")
  ))
}
#' ## test-drive_share.R
me_ <- nm_fun("TEST-drive-share")
nm_ <- nm_fun("TEST-drive-share", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("mirrors-to-share"),
    nm_("DESC")
  ))
}
#' ## test-drive_trash.R
me_ <- nm_fun("TEST-drive-trash")
nm_ <- nm_fun("TEST-drive-trash", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("trash-fodder"),
    me_("trashee-1"),
    me_("trashee-2")
  ))
}
#' ## test-drive_update.R
me_ <- nm_fun("TEST-drive-update")
nm_ <- nm_fun("TEST-drive-update", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("update-fodder"),
    nm_("not-unique"),
    nm_("does-not-exist")
  ))
}
#' ## test-drive_upload.R
me_ <- nm_fun("TEST-drive-upload")
nm_ <- nm_fun("TEST-drive-upload", NULL)
if (CLEAN) {
  drive_trash(c(
    nm_("upload-into-me"),
    nm_("DESCRIPTION")
  ))
}
