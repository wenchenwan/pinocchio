# ============================================================
# Pinocchio 源码架构三大机制 —— 可运行实验
#
# 用法：
#   make            编译全部三个实验
#   make run        依次运行全部实验
#   make run-01     只跑实验 1（三层文件布局）
#   make run-02     只跑实验 2（模板 + context 标量）
#   make run-03     只跑实验 3（CRTP 访问者）
#   make break-01   演示"删掉 .hxx include 会怎样"（实验 1 的反面教材）
#   make clean      清理
# ============================================================

CXX      ?= g++
CXXFLAGS ?= -std=c++14 -O2 -Wall

# 实验 2、3 需要链接 Pinocchio；实验 1 是自包含的玩具工程，不需要
PINO_FLAGS := $(shell pkg-config --cflags --libs pinocchio 2>/dev/null)

BUILD := build

.PHONY: all run run-01 run-02 run-03 break-01 clean check-deps

all: $(BUILD)/01-three-layer $(BUILD)/02-template-context $(BUILD)/03-visitor

check-deps:
	@if [ -z "$(PINO_FLAGS)" ]; then \
		echo "错误：pkg-config 找不到 pinocchio。"; \
		echo "请先安装 Pinocchio C++ 开发包，或设置 PKG_CONFIG_PATH。"; \
		exit 1; \
	fi

$(BUILD):
	@mkdir -p $(BUILD)

# ---- 实验 1：纯标准 C++，无需 Pinocchio ----
$(BUILD)/01-three-layer: 01-three-layer/main.cpp 01-three-layer/square.cpp \
                         01-three-layer/square.hpp 01-three-layer/square.hxx | $(BUILD)
	$(CXX) $(CXXFLAGS) -I01-three-layer \
		01-three-layer/main.cpp 01-three-layer/square.cpp -o $@

# ---- 实验 2、3：需要 Pinocchio ----
$(BUILD)/02-template-context: 02-template-context/main.cpp | $(BUILD)
	@$(MAKE) --no-print-directory check-deps
	$(CXX) $(CXXFLAGS) $< -o $@ $(PINO_FLAGS)

$(BUILD)/03-visitor: 03-visitor/main.cpp | $(BUILD)
	@$(MAKE) --no-print-directory check-deps
	$(CXX) $(CXXFLAGS) $< -o $@ $(PINO_FLAGS)

run: run-01 run-02 run-03

run-01: $(BUILD)/01-three-layer
	@./$<

run-02: $(BUILD)/02-template-context
	@./$<

run-03: $(BUILD)/03-visitor
	@./$<

# ---- 反面教材：删掉 .hpp 末尾的 include .hxx，观察链接失败 ----
break-01:
	@echo "=== 实验 1 反面教材：注释掉 square.hpp 末尾的 #include \"square.hxx\" ==="
	@rm -rf $(BUILD)/broken && mkdir -p $(BUILD)/broken
	@cp 01-three-layer/* $(BUILD)/broken/
	@sed -i 's|^#include "square.hxx"|// #include "square.hxx"  // ← 故意注释掉|' \
		$(BUILD)/broken/square.hpp
	@echo "--- 开始编译（预期：double 正常，float / long double 链接失败）---"
	@cd $(BUILD)/broken && $(CXX) $(CXXFLAGS) main.cpp square.cpp -o broken 2>&1 \
		| grep "undefined reference" | sed 's/.*undefined reference to /  链接失败: /' \
		|| echo "  （没有报错？请检查是否改对了文件）"
	@echo "--- 结论 ---"
	@echo "  double 版能链接，因为 square.cpp 已显式实例化并编进 .o"
	@echo "  float / long double 失败，因为没有 .hxx 就看不到模板定义体"
	@echo "  → 这就是 Pinocchio 每个 .hpp 末尾都要 include 对应 .hxx 的原因"

clean:
	@rm -rf $(BUILD)
	@echo "已清理 $(BUILD)/"
