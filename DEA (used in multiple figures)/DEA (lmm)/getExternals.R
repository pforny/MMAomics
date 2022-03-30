## to be extended to automatize as much as possible the retrieval  of external data.
if(!file.exists("csrproject")){
	system("git clone git@github.com:dlampart/csrproject.git")
	system("cd csrproject;git checkout smallEdits;cd ..")
}