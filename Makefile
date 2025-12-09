CXX := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -O2

SRC_DIR := src
BIN := app

# Все .cpp файлы в src/
SRCS := $(wildcard $(SRC_DIR)/*.cpp)

# Преобразуем список .cpp → .o
OBJS := $(SRCS:.cpp=.o)

all: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Правило компиляции каждого .cpp → .o
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Удаление временных файлов
clean:
	rm -f $(SRC_DIR)/*.o $(BIN)

# Полная очистка + бинарник
fclean: clean

re: fclean all

.PHONY: all clean fclean re
