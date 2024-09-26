Основная программа (анализ fastq) - main.py
Аргументы:

python3 main.py <input file> [--noask] [--export_fastq] [--export_csv] [--export_merge] [--export_dpi DPI]
    [--export FOLDER] [--adapters_len N] --mode [base, strip, plot] [--adapters FILE]

    input file - путь до fastq файла

    --noask - не спрашивать о длине адаптеров (сразу соглашаться)
    --export_fastq - экспортировать fastq файл, если были удалены адаптеры
    --export_csv - экспортировать статистику в файл

    --export_merge - отрисовывать графики одним файлом
    --export_dpi DPI - разрешение экспортируемых графиков
    --export FOLDER - папка, в которую будут экспортироваться файлы. Если не указано, ничего не экспортируется

    --adapters_len N - минимальная длина включения адаптеров. Может быть изменена интерактивно, если нет --noask
    --adapters FILE - файл с последовательностями адаптеров

    --mode - режим работы. Может быть несколько
        base - базовая статистика
        strip - удаление адаптеров. Если не указано, все аргументы адаптеров игнорируются
        plot - вывод графиков. Если не указано, все аргументы графиков игнорируются

Пример (для полного выполнения задания):
final_task_READS055722.student_7.fastq
--noask
--export_fastq
--export_csv
--export_merge
--adapters_len
5
--export_dpi
100
--export
./graphs/
--mode
base
strip
plot
--adapters
adapters.txt

Экспортируемые данные:
* fastq_sliced.fastq - fastq файл после отсечения последовательностей
* stats.csv - статистика:
Первые строки - property (параметр) и value (значение)
Далее - параметры отсечения:
    sequence (номер последовательности),
    offset_pos (отсечено с начала),
    offset_neg (отсечено с конца)

----

Программа-парсер GEO - main_geo.py
Аргументы:

python3 main_geo.py <URL>

Пример (из задания):
python3 main_geo.py "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM357351"