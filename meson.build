project('cpplibtest', 'cpp',
  default_options: ['cpp_std=c++20'])

# Указываем, что исходные файлы находятся в папке src
src_files = files('src/main.cpp')

# Создаём исполняемый файл
exe = executable('cpplibtest',
                 src_files,
                 install: false)  # Не устанавливаем в систему

configure_file(
  input: exe,
  output: 'cpplibtest',
  copy: true,
  install: false)