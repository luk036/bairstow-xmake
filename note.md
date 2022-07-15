# g++ vs clang

## structure binding

```cpp
auto [a, b, c] = get_result();
CHECK(a == 3);
```

## Range-based for loop

```cpp
for (auto a : foo.get_bar()) {
   // ...
}
```

## Concept

euqal_comparable
