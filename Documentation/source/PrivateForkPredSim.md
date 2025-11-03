Creating a private fork
=======================

The [PredSim repository](https://github.com/KULeuvenNeuromechanics/PredSim) is public and Github does not allow the creation of private forks for public repositories.
If you want to have a private version that can still fetch updates from the public repository, you can use this workaround.

Workflow below is adapted from [here](https://gist.github.com/0xjac/85097472043b697ab57ba1b1c7530274).

The correct way of creating a private frok by duplicating the repo is documented [here](https://help.github.com/articles/duplicating-a-repository/).

For this repository the Git Bash commands are:

 1. Create a bare clone of the repository.
    (This is temporary and will be removed so just do it wherever.)
    ```bash
    git clone --bare https://github.com/KULeuvenNeuromechanics/PredSim.git
    ```

 2. [Create a new private repository on Github](https://help.github.com/articles/creating-a-new-repository/) and name it `PredSim_private`.

 3. Mirror-push your bare clone to your new `PredSim_private` repository.
    > Replace `<your_username>` with your actual Github username in the url below.
    
    ```bash
    cd PredSim.git
    git push --mirror https://github.com/<your_username>/PredSim_private.git
    ```

 4. Remove the temporary local repository you created in step 1.
    ```bash
    cd ..
    rm -rf PredSim.git
    ```
    
 5. You can now clone your `PredSim_private` repository on your machine (for example in the `C:/GBW_MyPrograms` folder).
    ```bash
    cd C:/GBW_MyPrograms
    git clone https://github.com/<your_username>/PredSim_private.git
    ```
   
 6. If you want, add the original repo as remote to fetch (potential) future changes.
    ```bash
    cd PredSim_private
    git remote add upstream https://github.com/KULeuvenNeuromechanics/PredSim.git
    ```
    You can list all your remotes with `git remote -v`. You should see:
    ```
    origin	https://github.com/<your_username>/PredSim_private.git (fetch)
    origin	https://github.com/<your_username>/PredSim_private.git (push)
    upstream	https://github.com/KULeuvenNeuromechanics/PredSim_private.git (fetch)
    upstream	https://github.com/KULeuvenNeuromechanics/PredSim_private.git (push)
    ```
    > When you push, do so on `origin` with `git push origin`.
   
    > When you want to pull changes from `upstream` you can just fetch the remote and merge on top of your work.
    ```bash
    git fetch upstream
	git checkout local_branch_to_update
    git merge upstream/upstream_branch_to_update_from
    ```
    And solve the conflicts if any. You can also opt to do rebase instead of a merge, but this rewrites the git history. For more details on the difference between merge and rebase, see e.g.  https://stackoverflow.com/questions/16666089/whats-the-difference-between-Git-merge-and-git-rebase

 7. If you use GitHub Desktop, you have to add PredSim_private to the list of repositories:
   - Current repository
   - Add
   - Add existing repository...
   - Local path: `C:/GBW_MyPrograms/PredSim_private` (Use the same folder you selected in step 5.)
   - Add repository
   
 8. Your fork might have a different default branch. Set `master` or your own branch as default.

 