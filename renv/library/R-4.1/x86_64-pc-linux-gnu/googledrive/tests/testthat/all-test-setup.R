#' ---
#' title: googledrive test setup
#' date: '`r format(Sys.time())`'
#' output: html_document
#' ---

#' This script aggregates the test-related setup code from all test files.

library(googledrive)
source(testthat::test_path('helper.R'))
whoami <- drive_user()
whoami[c('displayName', 'emailAddress')]

## change this to TRUE when you are really ready to do this!
SETUP <- TRUE

#' ## test-drive_cp.R
me_ <- nm_fun("TEST-drive-cp")
nm_ <- nm_fun("TEST-drive-cp", NULL)
if (SETUP) {
  drive_mkdir(nm_("i-am-a-folder"))
  drive_mkdir(nm_("not-unique-folder"))
  drive_mkdir(nm_("not-unique-folder"))
  drive_upload(
    system.file("DESCRIPTION"),
    nm_("i-am-a-file")
  )
}
#' ## test-drive_create.R
me_ <- nm_fun("TEST-drive-create")
nm_ <- nm_fun("TEST-drive-create", NULL)
if (SETUP) {
  drive_mkdir(nm_("create-in-me"))
}
#' ## test-drive_download.R
nm_ <- nm_fun("TEST-drive-download", NULL)
if (SETUP) {
  drive_upload(system.file("DESCRIPTION"), nm_("DESC"))
  drive_upload(
    system.file("DESCRIPTION"),
    nm_("DESC-doc"),
    type = "document"
  )
}
#' ## test-drive_find.R
me_ <- nm_fun("TEST-drive-find")
nm_ <- nm_fun("TEST-drive-find", NULL)
if (SETUP) {
  drive_mkdir(nm_("find-me"))
  drive_upload(
    system.file("DESCRIPTION"),
    nm_("copy-me")
  )
}
#' ## test-drive_get.R
nm_ <- nm_fun("TEST-drive-get", NULL)
if (SETUP) {
  file_in_root <- drive_upload(
    system.file("DESCRIPTION"),
    name = nm_("thing01")
  )
  drive_upload(system.file("DESCRIPTION"), name = nm_("thing02"))
  drive_upload(system.file("DESCRIPTION"), name = nm_("thing03"))
  folder_in_root <- drive_mkdir(nm_("thing01"))
  folder_in_folder <- drive_mkdir(nm_("thing01"), parent = folder_in_root)
  file_in_folder_in_folder <- drive_cp(
    file_in_root,
    path = folder_in_folder,
    name = nm_("thing01")
  )
  drive_upload(
    system.file("DESCRIPTION"),
    path = folder_in_root,
    name = nm_("thing04")
  )

  folder_1_of_2 <- drive_mkdir(nm_("parent01"))
  folder_2_of_2 <- drive_mkdir(nm_("parent02"))
  child_of_2_parents <- drive_upload(
    system.file("DESCRIPTION"),
    path = folder_1_of_2,
    name = nm_("child_of_2_parents")
  )
  drive_update(child_of_2_parents, addParents = as_id(folder_2_of_2))
}
#' ## test-drive_ls.R
nm_ <- nm_fun("TEST-drive-ls", NULL)
if (SETUP) {
  drive_mkdir(nm_("list-me"))
  drive_upload(
    system.file("DESCRIPTION"),
    path = file.path(nm_("list-me"), nm_("DESCRIPTION"))
  )
  drive_upload(
    file.path(R.home("doc"), "html", "about.html"),
    path = file.path(nm_("list-me"), nm_("about-html"))
  )

  ## for testing `recursive = TRUE`
  top <- drive_mkdir(nm_("topdir"))
  drive_upload(
    system.file("DESCRIPTION"),
    path = top,
    name = nm_("apple"),
    type = "document",
    starred = TRUE
  )
  folder1_level1 <- drive_mkdir(nm_("folder1-level1"), parent = top)
  drive_mkdir(nm_("folder2-level1"), parent = top)
  x <- drive_upload(
    system.file("DESCRIPTION"),
    path = folder1_level1,
    name = nm_("banana"),
    type = "document"
  )
  folder1_level2 <- drive_mkdir(nm_("folder1-level2"), parent = folder1_level1)
  x <- drive_upload(
    system.file("DESCRIPTION"),
    path = folder1_level2,
    name = nm_("cranberry"),
    type = "document",
    starred = TRUE
  )
}
#' ## test-drive_mkdir.R
me_ <- nm_fun("TEST-drive-mkdir")
nm_ <- nm_fun("TEST-drive-mkdir", NULL)
if (SETUP) {
  drive_mkdir(nm_("OMNI-PARENT"))
}
#' ## test-drive_mv.R
me_ <- nm_fun("TEST-drive-mv")
nm_ <- nm_fun("TEST-drive-mv", NULL)
if (SETUP) {
  drive_mkdir(nm_("move-files-into-me"))
}
#' ## test-drive_publish.R
nm_ <- nm_fun("TEST-drive-publish", NULL)
if (SETUP) {
  drive_upload(
    file.path(R.home("doc"), "html", "about.html"),
    name = nm_("foo_doc"),
    type = "document"
  )
  drive_upload(
    file.path(R.home("doc"), "BioC_mirrors.csv"),
    name = nm_("foo_sheet"),
    type = "spreadsheet"
  )
  drive_upload(
    file.path(R.home("doc"), "html", "RLogo.pdf"),
    name = nm_("foo_pdf")
  )
}
#' ## test-drive_share.R
me_ <- nm_fun("TEST-drive-share")
nm_ <- nm_fun("TEST-drive-share", NULL)
if (SETUP) {
  drive_upload(system.file("DESCRIPTION"), nm_("DESC"))
}
#' ## test-drive_trash.R
me_ <- nm_fun("TEST-drive-trash")
nm_ <- nm_fun("TEST-drive-trash", NULL)
if (SETUP) {
  drive_upload(
    system.file("DESCRIPTION"),
    nm_("trash-fodder")
  )
}
#' ## test-drive_update.R
me_ <- nm_fun("TEST-drive-update")
nm_ <- nm_fun("TEST-drive-update", NULL)
if (SETUP) {
  drive_upload(system.file("DESCRIPTION"), nm_("update-fodder"))
  drive_upload(system.file("DESCRIPTION"), nm_("not-unique"))
  drive_upload(system.file("DESCRIPTION"), nm_("not-unique"))
}
#' ## test-drive_upload.R
me_ <- nm_fun("TEST-drive-upload")
nm_ <- nm_fun("TEST-drive-upload", NULL)
if (SETUP) {
  drive_mkdir(nm_("upload-into-me"))
  drive_mkdir(nm_("upload-into-me-too"))
}
