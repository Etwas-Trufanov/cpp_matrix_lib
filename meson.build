project('cpplibtest', 'cpp',
  default_options: ['cpp_std=c++23'])

# Указываем, что исходные файлы находятся в папке src
src_files = files('src/main.cpp')

# Создаём исполняемый файл
exe = executable('cpplibtest',
                 src_files,
                 install: false)  # Не устанавливаем в систему
