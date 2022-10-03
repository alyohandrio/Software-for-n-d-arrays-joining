# Software for n-d arrays joining

Этот репозиторий содержит библиотеку, предназначенную для объединения многомерных неиндексированных массивов, хранящих значения типа `float`.


Класс, используемый для хранения массивов называется `Database`. Для инициализации его данными используется метод `simlend`, принимающий имя файла с данными, вектор размерностей массива, вектор предлагаемых размеров чанков, наименьшее значение, наибольшее значение, количество бинов, на которые разбивается гиперпрямоугольная область, соответствующая одному чанку, количество бинов в одном чанке (это число должно делить предыдущий параметр).

Для соединения массивов по признаку совпадения координат и соответствующих значений используется функция `equal_join`, принимающая два массива. Параметры, переданные в `simplend` для этих массивов должны совпадать.

Для соединения массивов по признаку совпадения координат и отличия соответсвующих значений на меньшее определённой величины число испольуется функция `value_similarity_join`. Функция принимает два массива, функцию сравнения и допустимую разницу. Ожидается, что функция сравнения нестрого возрастает с ростом модуля разности значений.