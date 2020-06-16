# TRD-self-tracking

Please remember to use `git pull`, before `git push` on order to avoid merge commits.

### Coding style

We can use the same style file as the ALICE O2 project. It is added to the repository now and in order to format your code accordingly you need to run
> clang-format -style=file YOUR_FILE

See the readme from the project here: https://github.com/AliceO2Group/CodingGuidelines

`clang-format` is installed on the Heidelberg servers. In order to use it you need to load clang via
> alienv load VO_ALICE@Clang::latest
