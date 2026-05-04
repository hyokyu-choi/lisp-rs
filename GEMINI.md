# Vibe Coding Workflow

This project follows a strict **Command-Driven TDD (Test-Driven Development)** workflow. All agents must adhere to these mandates.

## 1. Execution Cycle
When a directive (command) is received, follow this exact sequence:
1.  **Test First**: Create or update a test case that defines the expected behavior or reproduces the reported bug.
2.  **Implement/Fix**: Modify the codebase to make the test pass.
3.  **Validate**: Run the tests.
4.  **Iterate**: If tests fail, diagnose the issue, modify the code, and re-run tests until they pass.
5.  **Final Verification**: Run project-wide checks (lint, types) if applicable.

## 2. Strict Scope Control
- **No Unsolicited Actions**: Do not perform refactorings, cleanups, or optimizations outside the immediate scope of the command.
- **No Proactive Suggestions**: Do not offer recommendations, alternative approaches, or "just-in-case" code unless explicitly asked for an inquiry or advice.
- **Minimal Guessing**: If a command is underspecified, ask for clarification instead of making assumptions that lead to unintended code generation.
- **No Postambles**: Once a task is verified and complete, provide a concise confirmation. Avoid conversational filler or summaries of future possibilities.

## 3. Engineering Standards
- Maintain existing architectural patterns (Trait-based systems in `math/integrate.rs`, etc.).
- Ensure all new features are accompanied by unit tests in the same file or a corresponding `tests/` file.
- Use parallel tool execution whenever safe to minimize turn count.
