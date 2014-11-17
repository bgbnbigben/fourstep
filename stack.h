// Copyright 2011 - 2014 Brian Marshall. All rights reserved.
//
// Use of this source code is governed by the BSD 2-Clause License that can be
// found in the LICENSE file.

#ifndef SHUNTING_YARD_STACK_H
#define SHUNTING_YARD_STACK_H

typedef struct Stack Stack;

// Inserts a new element at the top of the stack. Initialize new stacks to NULL
// before calling this function.
void stack_push(Stack **stack, const void *value);

// Removes an element from the top of the stack. If empty, returns NULL.
const void *stack_pop(Stack **stack);

// Returns the value at the top of the stack.
const void *stack_top(const Stack *stack);

#endif  // SHUNTING_YARD_STACK_H
